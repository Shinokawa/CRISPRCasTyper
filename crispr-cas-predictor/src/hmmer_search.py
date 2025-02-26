import os
import subprocess
import tempfile
import logging
import multiprocessing
import sys  # 添加sys模块导入
import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO

class HMMERSearch:
    def __init__(self, hmm_model_dir, num_threads=None):
        """初始化HMMER搜索类
        
        Args:
            hmm_model_dir: 包含HMM模型的目录路径
            num_threads: 并行处理的线程数，默认使用可用CPU核心数
        """
        self.hmm_model_dir = hmm_model_dir
        if not os.path.exists(hmm_model_dir):
            raise FileNotFoundError(f"HMM模型目录不存在: {hmm_model_dir}")
            
        self.hmm_files = [os.path.join(hmm_model_dir, f) for f in os.listdir(hmm_model_dir) 
                         if f.endswith('.hmm')]
        
        if not self.hmm_files:
            raise FileNotFoundError(f"在 {hmm_model_dir} 中未找到HMM模型文件")
        
        # 设置线程数
        self.num_threads = num_threads if num_threads else multiprocessing.cpu_count()
        # 每个HMMER进程使用的CPU数量
        self.hmmer_cpu = max(1, min(2, self.num_threads // 4))
            
        logging.info(f"找到 {len(self.hmm_files)} 个HMM模型文件，使用 {self.num_threads} 个处理线程")
        
    def search(self, fasta_file):
        """搜索FASTA文件中的序列
        
        Args:
            fasta_file: 输入的FASTA文件路径
        
        Returns:
            包含HMMER搜索结果的列表
        """
        if not os.path.exists(fasta_file):
            raise FileNotFoundError(f"FASTA文件不存在: {fasta_file}")
            
        results = []
        logging.info(f"开始使用 {len(self.hmm_files)} 个HMM模型并行搜索序列...")
        
        # 计算并行搜索的最大工作进程数
        max_workers = min(self.num_threads, len(self.hmm_files))
        
        # 使用并行处理执行HMMER搜索
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # 提交所有任务
            future_to_hmm = {
                executor.submit(self.run_hmmer, fasta_file, hmm_file): 
                os.path.basename(hmm_file).replace('.hmm', '') 
                for hmm_file in self.hmm_files
            }
            
            # 使用tqdm创建进度条
            total = len(self.hmm_files)
            with tqdm.tqdm(total=total, desc="HMM搜索进度") as pbar:
                # 处理完成的任务
                for future in as_completed(future_to_hmm):
                    hmm_name = future_to_hmm[future]
                    try:
                        hmm_results = future.result()
                        if hmm_results:
                            parsed_results = self.parse_results(hmm_results, hmm_name)
                            results.extend(parsed_results)
                    except Exception as exc:
                        logging.error(f"处理 {hmm_name} 时出错: {exc}")
                    finally:
                        pbar.update(1)
                
        logging.info(f"HMMER搜索完成，找到 {len(results)} 个匹配")
        return results
    
    def search_chunks(self, chunk_files):
        """搜索多个FASTA文件块
        
        Args:
            chunk_files: FASTA文件块路径列表
        
        Returns:
            包含HMMER搜索结果的列表
        """
        all_results = []
        
        for i, chunk_file in enumerate(chunk_files):
            logging.info(f"处理FASTA块 {i+1}/{len(chunk_files)}: {os.path.basename(chunk_file)}")
            chunk_results = self.search(chunk_file)
            all_results.extend(chunk_results)
            
            # 搜索完成后删除临时块文件（如果文件名包含.chunk）
            if '.chunk' in chunk_file:
                try:
                    os.unlink(chunk_file)
                    logging.debug(f"已删除临时块文件: {chunk_file}")
                except Exception as e:
                    logging.warning(f"删除临时块文件时出错: {e}")
        
        return all_results

    def run_hmmer(self, fasta_file, hmm_file):
        """运行HMMER搜索
        
        Args:
            fasta_file: 输入的FASTA文件路径
            hmm_file: HMM模型文件路径
        
        Returns:
            HMMER命令的输出
        """
        try:
            with tempfile.NamedTemporaryFile(suffix='.out', delete=False) as temp_output:
                temp_output_path = temp_output.name
                
            hmm_name = os.path.basename(hmm_file)
            
            if sys.platform == "win32":
                # Windows平台使用Docker
                abs_fasta = os.path.abspath(fasta_file)
                abs_hmm = os.path.abspath(hmm_file)
                abs_output = os.path.abspath(temp_output_path)
                
                fasta_dir = os.path.dirname(abs_fasta)
                hmm_dir = os.path.dirname(abs_hmm)
                output_dir = os.path.dirname(abs_output)
                
                cmd = ['docker', 'run', '--rm',
                      f'-v', f'{fasta_dir}:/fasta',
                      f'-v', f'{hmm_dir}:/hmm',
                      f'-v', f'{output_dir}:/output',
                      'hmmer-image',
                      'hmmsearch',
                      '--cpu', str(self.hmmer_cpu),
                      '--tblout', f'/output/{os.path.basename(temp_output_path)}',
                      '--noali',  # 不输出比对结果以节省内存
                      '-Z', '50000000',  # 调整数据库大小假设值
                      '--F1', '0.02',  # MSV过滤器阈值
                      '--F2', '0.001',  # Viterbi过滤器阈值
                      f'/hmm/{os.path.basename(hmm_file)}',
                      f'/fasta/{os.path.basename(fasta_file)}']
            else:
                # Unix平台直接运行HMMER
                cmd = ['hmmsearch', 
                      '--cpu', str(self.hmmer_cpu), 
                      '--tblout', temp_output_path,
                      '--noali',  # 不输出比对结果以节省内存
                      '-Z', '50000000',  # 调整数据库大小假设值
                      '--F1', '0.02',  # MSV过滤器阈值
                      '--F2', '0.001',  # Viterbi过滤器阈值
                      hmm_file, 
                      fasta_file]
            
            process = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            if os.path.exists(temp_output_path):
                with open(temp_output_path, 'r') as f:
                    lines = f.readlines()
                os.unlink(temp_output_path)  # 删除临时文件
                return lines
            return []
            
        except subprocess.CalledProcessError as e:
            logging.error(f"HMMER搜索失败: {e}")
            logging.error(f"标准输出: {e.stdout.decode('utf-8') if hasattr(e, 'stdout') else ''}")
            logging.error(f"标准错误: {e.stderr.decode('utf-8') if hasattr(e, 'stderr') else ''}")
            return []
        except Exception as e:
            logging.error(f"运行HMMER时出错: {e}")
            return []

    def parse_results(self, hmmer_output, hmm_name):
        """解析HMMER输出
        
        Args:
            hmmer_output: HMMER命令的输出行
            hmm_name: HMM模型的名称
        
        Returns:
            解析后的结果列表
        """
        results = []
        
        for line in hmmer_output:
            if line.startswith('#') or not line.strip():
                continue
                
            try:
                parts = line.split(None, 22)  # 使用空白分隔，最多分成23部分
                
                if len(parts) < 9:  # 最低需要9个字段
                    continue
                    
                # tblout格式: target_name accession query_name accession e-value score bias ...
                query_name = parts[0]  # 蛋白质序列ID
                target_name = parts[2]  # HMM模型名称
                
                # 如果没有提供具体的目标名称，使用HMM文件名
                if target_name == '-':
                    target_name = hmm_name
                
                try:
                    e_value = float(parts[4])
                except ValueError:
                    e_value = float('inf')
                    
                try:
                    score = float(parts[5])
                except ValueError:
                    score = 0.0
                
                result = {
                    'query_name': query_name,
                    'target_name': target_name,
                    'e_value': e_value,
                    'score': score,
                    'hmm_model': hmm_name
                }
                
                results.append(result)
                
            except Exception as e:
                logging.warning(f"解析HMMER输出行时出错: {line.strip()}")
                logging.warning(f"错误: {e}")
                
        return results