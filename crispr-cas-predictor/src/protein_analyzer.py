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
        try:
            self.sequences = list(SeqIO.parse(input_file, "fasta"))
            logging.info(f"从{input_file}加载了{len(self.sequences)}个序列")
        except Exception as e:
            logging.error(f"无法解析FASTA文件: {e}")
            self.sequences = []

    def filter_sequences(self, min_length=50):
        """根据最小长度过滤序列
        
        Args:
            min_length: 最小蛋白质长度
            
        Returns:
            过滤后的序列列表
        """
        original_count = len(self.sequences)
        self.sequences = [seq for seq in self.sequences if len(seq.seq) >= min_length]
        filtered_count = original_count - len(self.sequences)
        logging.info(f"过滤掉了{filtered_count}个短于{min_length}氨基酸的序列")
        return self.sequences

    def preprocess_sequences(self):
        """预处理序列（例如，移除重复项，标准化）
        
        Returns:
            预处理后的序列列表
        """
        # 在这里可以添加更多的预处理步骤
        seen_ids = set()
        unique_sequences = []
        
        for seq in self.sequences:
            if seq.id not in seen_ids:
                seen_ids.add(seq.id)
                unique_sequences.append(seq)
        
        duplicates_count = len(self.sequences) - len(unique_sequences)
        if duplicates_count > 0:
            logging.info(f"移除了{duplicates_count}个重复序列")
            self.sequences = unique_sequences
            
        return self.sequences

    def write_filtered_sequences(self, output_file=None):
        """将过滤后的序列写入文件
        
        Args:
            output_file: 输出文件路径，如果为None则创建临时文件
            
        Returns:
            包含过滤序列的文件路径
        """
        if output_file is None:
            with tempfile.NamedTemporaryFile(suffix='.fasta', delete=False) as temp_file:
                output_file = temp_file.name
                
        SeqIO.write(self.sequences, output_file, "fasta")
        logging.info(f"将{len(self.sequences)}个序列写入{output_file}")
        return output_file

    def analyze(self):
        """对蛋白质序列进行分析
        
        Returns:
            分析结果
        """
        self.filter_sequences()
        self.preprocess_sequences()
        
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