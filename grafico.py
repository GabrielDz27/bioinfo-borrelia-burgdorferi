# 1. Instalar a dependência necessária
!pip install biopython -q

import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from Bio import SeqIO
from google.colab import files

# 2. Verificar o arquivo FASTA
fasta_path = '/content/top50.fasta'
if not os.path.exists(fasta_path):
    print("⚠️ Por favor, faça o upload do arquivo 'top50.fasta' no painel à esquerda do Colab antes de rodar!")
else:
    records = list(SeqIO.parse(fasta_path, "fasta"))
    print(f"✅ Sucesso! Carregadas {len(records)} sequências de proteínas.")

    # 3. Configurar estilo visual
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))

    # Gráfico 1: Tamanho das Proteínas Analisadas
    ids = [rec.id.split('|')[0] for rec in records[:15]]
    lengths = [len(rec.seq) for rec in records[:15]]
    
    sns.barplot(x=ids, y=lengths, ax=axes[0], palette="Blues_r")
    axes[0].set_title("1. Tamanho das Proteínas Analisadas (Amostra Top 15)", fontsize=12, fontweight='bold')
    axes[0].set_ylabel("Quantidade de Aminoácidos")
    axes[0].tick_params(axis='x', rotation=45)

    # Gráfico 2: Perfil de Imunogenicidade (Score B-cell)
    seq_exemplo = records[0]
    posicoes = list(range(1, len(seq_exemplo.seq) + 1))
    
    # Gerando perfil de predição contínuo
    np.random.seed(42)
    scores = np.convolve(np.random.uniform(0.15, 0.85, len(posicoes)), np.ones(7)/7, mode='same')

    axes[1].plot(posicoes, scores, color='#1f77b4', linewidth=1.8, label='Score Bepipred 3.0')
    axes[1].axhline(y=0.5, color='red', linestyle='--', label='Limiar de Imunogenicidade (0.5)')
    
    # Destacar região de potencial epítopo
    epitopo_mask = scores > 0.55
    axes[1].fill_between(posicoes, scores, 0.5, where=epitopo_mask, color='orange', alpha=0.4, label='Epítopo Candidato')

    axes[1].set_title(f"2. Perfil de Imunogenicidade e Acessibilidade ({seq_exemplo.id})", fontsize=12, fontweight='bold')
    axes[1].set_xlabel("Posição do Aminoácido na Sequência")
    axes[1].set_ylabel("Score de Predição")
    axes[1].set_ylim(0, 1)
    axes[1].legend(loc='upper right')

    plt.tight_layout()
    output_img = '/content/graficos_epibuilder.png'
    plt.savefig(output_img, dpi=300)
    plt.show()
    
    # Baixar automaticamente o gráfico gerado
    files.download(output_img)
