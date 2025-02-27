from Bio import SeqIO
import os
import tempfile
import logging
import gzip

class ProteinAnalyzer:
    def __init__(self, input_file, estimate_count=True):
        """初始化蛋白质分析器
        
        Args:
            input_file: 输入的FASTA文件路径
            estimate_count: 是否估计序列数量（对于大文件设为False可提高性能）
        """
        self.input_file = input_file
        self.file_size = os.path.getsize(input_file) if os.path.exists(input_file) else 0
        
        # 对于大文件，使用文件大小估算序列数量
        if self.file_size > 1024 * 1024 * 1000:  # 超过1GB
            # 假设每个序列平均200个氨基酸，每个字符1字节，再加上FASTA头部约50字节
            estimated_seq_size = 250  # 字节
            self.sequence_count = self.file_size // estimated_seq_size
            logging.info(f"大型FASTA文件检测到 ({self.file_size/(1024*1024*1024):.2f} GB)，估计包含约 {self.sequence_count} 个序列")
        elif estimate_count:
            self.sequence_count = self._estimate_sequences()
            logging.info(f"FASTA文件包含约 {self.sequence_count} 个序列")
        else:
            self.sequence_count = None
    
    def _estimate_sequences(self):
        """估计FASTA文件中的序列数量
        
        对于大文件，只采样前100MB来估计
        """
        try:
            # 检查是否为压缩文件
            is_gzip = self.input_file.endswith('.gz')
            
            count = 0
            bytes_read = 0
            max_sample = 100 * 1024 * 1024  # 100MB采样
            
            if is_gzip:
                opener = gzip.open
            else:
                opener = open
                
            with opener(self.input_file, 'rt') as f:
                for line in f:
                    bytes_read += len(line)
                    if line.startswith('>'):
                        count += 1
                    if bytes_read > max_sample:
                        break
            
            # 如果只读了部分文件，根据比例估计总数
            if bytes_read < self.file_size and bytes_read > 0:
                count = int(count * (self.file_size / bytes_read))
                
            return count
        except Exception as e:
            logging.error(f"估计序列数量时出错: {e}")
            return 0
    
    def filter_sequences(self, min_length=50, output_file=None):
        """流式过滤序列并直接写入输出文件
        
        Args:
            min_length: 最小蛋白质长度
            output_file: 输出文件路径，如果为None则创建临时文件
            
        Returns:
            包含过滤序列的文件路径
        """
        if output_file is None:
            with tempfile.NamedTemporaryFile(suffix='.fasta', delete=False) as temp_file:
                output_file = temp_file.name
                
        filtered_count = 0
        total_count = 0
        duplicate_ids = set()
        
        # 检查是否为压缩文件
        is_gzip = self.input_file.endswith('.gz')
        
        # 使用流处理方式读取和写入
        try:
            # 使用适当的打开方式
            if is_gzip:
                input_handle = gzip.open(self.input_file, "rt")
            else:
                input_handle = open(self.input_file, "r")
                
            with input_handle, open(output_file, 'w') as out_f:
                # 批量处理序列
                batch_size = 5000  # 调整批次大小以平衡内存使用和性能
                current_batch = []
                current_batch_size = 0
                
                for record in SeqIO.parse(input_handle, "fasta"):
                    total_count += 1
                    
                    # 每处理10000条序列输出一次进度
                    if total_count % 10000 == 0:
                        logging.info(f"已处理 {total_count} 个序列，过滤了 {filtered_count} 个序列")
                    
                    # 长度过滤
                    if len(record.seq) < min_length:
                        filtered_count += 1
                        continue
                        
                    # 重复ID过滤（可选，对大型数据集可能需要禁用）
                    if record.id in duplicate_ids:
                        filtered_count += 1
                        continue
                    
                    duplicate_ids.add(record.id)
                    
                    # 添加到当前批次
                    current_batch.append(record)
                    current_batch_size += 1
                    
                    # 批次满了就写入文件并清空批次
                    if current_batch_size >= batch_size:
                        SeqIO.write(current_batch, out_f, "fasta")
                        current_batch = []
                        current_batch_size = 0
                        
                    # 定期清空重复ID集合以控制内存使用
                    if len(duplicate_ids) > 100000:
                        duplicate_ids.clear()  # 完全清空以节省内存
                
                # 写入最后一批
                if current_batch:
                    SeqIO.write(current_batch, out_f, "fasta")
            
        except Exception as e:
            logging.error(f"过滤序列时出错: {e}")
            import traceback
            logging.error(traceback.format_exc())
            raise
        
        logging.info(f"过滤完成: 总共处理 {total_count} 个序列，过滤掉 {filtered_count} 个序列")
        logging.info(f"将过滤后的序列写入 {output_file}")
        return output_file
    
    def preprocess_and_write(self, output_file=None, min_length=50, batch_size=1000):
        """流式预处理序列并写入文件
        
        统一所有处理步骤成一个流程，只读一次文件
        
        Args:
            output_file: 输出文件路径，如果为None则创建临时文件
            min_length: 最小蛋白质长度
            batch_size: 批处理大小
            
        Returns:
            包含预处理序列的文件路径
        """
        return self.filter_sequences(min_length, output_file)
    
    def analyze(self):
        """对蛋白质序列进行分析
        
        Returns:
            分析结果
        """
        self.filter_sequences()
        
        # 基本序列统计
        stats = {
            "total_sequences": len(self.sequences),
            "avg_length": sum(len(seq.seq) for seq in self.sequences) / len(self.sequences) if self.sequences else 0,
            "min_length": min(len(seq.seq) for seq in self.sequences) if self.sequences else 0,
            "max_length": max(len(seq.seq) for seq in self.sequences) if self.sequences else 0
        }
        
        logging.info(f"序列分析完成: {stats['total_sequences']}个序列, "
                    f"平均长度: {stats['avg_length']:.2f}, "
                    f"最短: {stats['min_length']}, "
                    f"最长: {stats['max_length']}")
                    
        return stats

    def split_fasta_file(self, chunk_size=1000000000):  # 约1GB每块
        """将大型FASTA文件分割成多个小文件
        
        Args:
            chunk_size: 每个文件的近似大小(字节)
                
        Returns:
            分割后的文件路径列表
        """
        chunk_files = []
        
        try:
            # 如果文件小于分块大小，不做分块
            if self.file_size <= chunk_size:
                logging.info(f"文件大小 ({self.file_size/(1024*1024):.2f} MB) 小于分块大小，不需要分割")
                return [self.input_file]
                
            logging.info(f"文件大小为 {self.file_size/(1024*1024*1024):.2f} GB，进行分块处理")
            
            # 检查是否为压缩文件
            is_gzip = self.input_file.endswith('.gz')
            
            # Windows平台特定优化: 直接按行读取，避免使用SeqIO解析全部记录
            chunk_index = 0
            current_size = 0
            buffer_size = 10 * 1024 * 1024  # 10MB写入缓冲区
            
            # 创建第一个输出文件
            chunk_file = f"{self.input_file}.chunk{chunk_index}.fasta"
            out_f = open(chunk_file, 'w', buffering=buffer_size)
            chunk_files.append(chunk_file)
            
            # 打开输入文件
            if is_gzip:
                input_handle = gzip.open(self.input_file, "rt")
            else:
                input_handle = open(self.input_file, "r")
                
            with input_handle:
                sequence = ""
                header = ""
                for line in input_handle:
                    line_size = len(line)
                    
                    if line.startswith('>'):
                        # 处理之前的序列
                        if header and sequence:
                            seq_size = len(header) + len(sequence) + 2  # 加上换行符
                            
                            # 如果当前块已满，创建新块
                            if current_size > 0 and current_size + seq_size > chunk_size:
                                out_f.close()
                                logging.info(f"写入分块文件 {chunk_file} 完成 ({current_size/(1024*1024):.2f} MB)")
                                
                                # 创建新文件
                                chunk_index += 1
                                chunk_file = f"{self.input_file}.chunk{chunk_index}.fasta"
                                out_f = open(chunk_file, 'w', buffering=buffer_size)
                                chunk_files.append(chunk_file)
                                current_size = 0
                            
                            # 写入序列到当前文件
                            out_f.write(f"{header}\n{sequence}\n")
                            current_size += seq_size
                        
                        # 保存新的序列头
                        header = line.strip()
                        sequence = ""
                    else:
                        # 累积序列
                        sequence += line.strip()
                
                # 处理最后一个序列
                if header and sequence:
                    out_f.write(f"{header}\n{sequence}\n")
            
            # 关闭最后一个文件
            out_f.close()
            logging.info(f"写入最后一个分块文件 {chunk_file} 完成")
            
            return chunk_files
            
        except Exception as e:
            logging.error(f"分块处理文件时出错: {e}")
            import traceback
            logging.error(traceback.format_exc())
            
            # 确保所有打开的文件都被关闭
            try:
                if 'out_f' in locals() and not out_f.closed:
                    out_f.close()
            except:
                pass
                
            return [self.input_file]  # 出错时返回原文件