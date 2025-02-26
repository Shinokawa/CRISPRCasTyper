import os
import json
import pandas as pd
import logging

class CASTyping:
    def __init__(self, identified_proteins, cas_types_file=None):
        """初始化CAS分型类
        
        Args:
            identified_proteins: HMMER搜索结果
            cas_types_file: 包含蛋白质类型信息的文件（可选）
        """
        self.identified_proteins = identified_proteins
        self.cas_types_file = cas_types_file
        self.cas_type_mapping = self._load_cas_type_mapping()
        
    def _load_cas_type_mapping(self):
        """加载Cas蛋白类型映射
        
        从文件或内部定义中加载Cas蛋白类型与HMM名称的映射关系
        """
        # 默认映射，基于hmm文件名的关键字
        default_mapping = {
            # Type I
            'cas3': 'Type I (Cas3)',
            'cas8': 'Type I',
            'cas5': 'Type I',
            'cas6': 'Type I',
            'cas7': 'Type I',
            'cse1': 'Type I-E',
            'cse2': 'Type I-E',
            'csy1': 'Type I-F',
            'csy2': 'Type I-F',
            'csy3': 'Type I-F',
            'csa5': 'Type I-A',
            'csm': 'Type III-A',
            # Type II
            'cas9': 'Type II (Cas9)',
            'csn2': 'Type II-A',
            # Type III
            'cas10': 'Type III (Cas10)',
            'cmr': 'Type III-B',
            'csx': 'Type III',
            # Type IV
            'csf': 'Type IV',
            # Type V
            'cas12': 'Type V (Cas12)',
            'cas12a': 'Type V-A',
            'cas12b': 'Type V-B',
            'cas12c': 'Type V-C',
            # Type VI
            'cas13': 'Type VI (Cas13)',
            'cas13a': 'Type VI-A',
            'cas13b': 'Type VI-B',
            'cas13c': 'Type VI-C',
            'cas13d': 'Type VI-D',
            # Common Cas proteins
            'cas1': 'Adaptation',
            'cas2': 'Adaptation',
            'cas4': 'Adaptation',
            'dinG': 'Adaptation',
        }
        
        # 如果提供了专门的映射文件，则从文件加载
        if self.cas_types_file and os.path.exists(self.cas_types_file):
            try:
                with open(self.cas_types_file, 'r') as f:
                    mapping = json.load(f)
                return mapping
            except Exception as e:
                logging.warning(f"无法从文件加载CAS类型映射: {e}，使用默认映射")
                
        return default_mapping

    def classify_proteins(self):
        """分类识别到的蛋白质
        
        Returns:
            一个字典，将蛋白质ID映射到其预测的类型
        """
        classified_results = {}
        
        for protein in self.identified_proteins:
            protein_id = protein['query_name']
            hmm_model = protein.get('target_name', '').lower()
            e_value = protein['e_value']
            score = protein.get('score', 0)
            
            # 确定置信度
            if e_value <= 1e-10:
                confidence = "High"
            elif e_value <= 1e-5:
                confidence = "Medium"
            else:
                confidence = "Low"
            
            # 计算概率得分 (0-100区间)
            # 使用E-value的负对数来表示概率
            # 较低的E-value会得到较高的概率分数
            import math
            if e_value > 0:
                prob_score = min(100, -10 * math.log10(e_value))
            else:
                prob_score = 100
            
            # 确定CRISPR-Cas类型
            cas_type = "Unknown"
            for key, type_value in self.cas_type_mapping.items():
                if key.lower() in hmm_model.lower():
                    cas_type = type_value
                    break
            
            # 创建结果信息，包含概率得分
            result = {
                'type': cas_type,
                'confidence': confidence,
                'e_value': e_value,
                'score': score,
                'probability': prob_score,
                'hmm_model': hmm_model
            }
            
            classified_results[protein_id] = result
            
        return classified_results

    def run_typing(self):
        """运行完整的分型流程"""
        return self.classify_proteins()
    
    def generate_report(self, output_file=None):
        """生成分型报告
        
        Args:
            output_file: 输出文件路径
        """
        results = self.classify_proteins()
        
        # 统计每种类型的数量
        type_counts = {}
        for protein_id, result in results.items():
            cas_type = result['type']
            if cas_type not in type_counts:
                type_counts[cas_type] = 0
            type_counts[cas_type] += 1
        
        # 创建报告
        report = ["# CRISPR-Cas Protein Classification Report", ""]
        report.append("## Summary")
        report.append("Number of proteins analyzed: {}".format(len(results)))
        report.append("")
        
        report.append("## Type Distribution")
        for cas_type, count in type_counts.items():
            report.append("- {}: {} proteins".format(cas_type, count))
        report.append("")
        
        # 按概率分数排序结果
        sorted_results = sorted(
            [(pid, res) for pid, res in results.items()], 
            key=lambda x: x[1]['probability'], 
            reverse=True
        )
        
        report.append("## Detailed Results (Sorted by Probability)")
        for protein_id, result in sorted_results:
            report.append("- {}: {} (Confidence: {}, Probability: {:.1f}%, E-value: {:.4g})".format(
                protein_id, 
                result['type'], 
                result['confidence'], 
                result['probability'],
                result['e_value']
            ))
        
        # 保存或打印报告
        if output_file:
            with open(output_file, 'w') as f:
                f.write('\n'.join(report))
        
        return '\n.join(report)'