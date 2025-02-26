import os
import tempfile
from Bio import SeqIO
import logging

class ProteinAnalyzer:
    def __init__(self, input_file):
        """初始化蛋白质分析器
        
        Args:
            input_file: 输入的FASTA文件路径
        """
        self.input_file = input_file
        self.sequence_count = self._count_sequences()
        logging.info(f"FASTA文件包含约 {self.sequence_count} 个序列")
    
    def _count_sequences(self):
        """计算FASTA文件中的序列数量"""
        count = 0
        try:
            with open(self.input_file) as f:
                for line in f:
                    if line.startswith('>'):
                        count += 1
            return count
        except Exception as e:
            logging.error(f"计算序列数量时出错: {e}")
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
        
        # 使用流处理方式读取和写入
        with open(output_file, 'w') as out_f:
            # 批量处理序列
            batch_size = 1000  # 调整批次大小以平衡内存使用和性能
            current_batch = []
            
            for record in SeqIO.parse(self.input_file, "fasta"):
                total_count += 1
                
                # 长度过滤
                if len(record.seq) < min_length:
                    filtered_count += 1
                    continue
                    
                # 重复ID过滤
                if record.id in duplicate_ids:
                    continue
                duplicate_ids.add(record.id)
                
                # 添加到当前批次
                current_batch.append(record)
                
                # 批次满了就写入文件并清空批次
                if len(current_batch) >= batch_size:
                    SeqIO.write(current_batch, out_f, "fasta")
                    current_batch = []
                    
                # 定期清空重复ID集合以控制内存使用
                if len(duplicate_ids) > 100000:  # 调整为适当的值
                    duplicate_ids = set()
            
            # 写入最后一批
            if current_batch:
                SeqIO.write(current_batch, out_f, "fasta")
        
        logging.info(f"过滤掉了 {filtered_count} 个短于 {min_length} 氨基酸的序列，共处理 {total_count} 个序列")
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
            file_size = os.path.getsize(self.input_file)
            if file_size <= chunk_size:
                return [self.input_file]  # 文件小于分块大小，不需要分割
                
            logging.info(f"文件大小为 {file_size/(1024*1024*1024):.2f} GB，进行分块处理")
            
            chunk_index = 0
            current_size = 0
            records = []
            
            for record in SeqIO.parse(self.input_file, "fasta"):
                # 估算记录大小（粗略近似）
                record_size = len(str(record.seq)) + len(record.id) + 50
                
                # 如果当前块将超出最大大小，则写入文件
                if current_size + record_size > chunk_size and records:
                    chunk_file = f"{self.input_file}.chunk{chunk_index}.fasta"
                    with open(chunk_file, 'w') as out_f:
                        SeqIO.write(records, out_f, "fasta")
                    chunk_files.append(chunk_file)
                    records = []
                    current_size = 0
                    chunk_index += 1
                    logging.info(f"写入分块文件 {chunk_file}")
                
                records.append(record)
                current_size += record_size
            
            # 写入最后一个块
            if records:
                chunk_file = f"{self.input_file}.chunk{chunk_index}.fasta"
                with open(chunk_file, 'w') as out_f:
                    SeqIO.write(records, out_f, "fasta")
                chunk_files.append(chunk_file)
                logging.info(f"写入最后一个分块文件 {chunk_file}")
            
            return chunk_files
            
        except Exception as e:
            logging.error(f"分块处理文件时出错: {e}")
            return [self.input_file]  # 出错时返回原文件