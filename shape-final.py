import pandas as pd
import numpy as np
import argparse

def process_promoter(promoter, tss_positions):
    seqname = promoter['seqnames']
    start = promoter['start']
    end = promoter['end']
    
    print(f"Processing promoter: {seqname} {start} {end}")  
    
    filtered_tss = tss_positions[(tss_positions['seqnames'] == seqname) & 
                                 (tss_positions['pos'] >= start) & 
                                 (tss_positions['pos'] <= end)]
    
    if filtered_tss.empty:
        return pd.DataFrame({
            'seqnames': [seqname],
            'start': [start],
            'end': [end],
            'width': [end - start + 1],
            'total_tss': [0],
            'tss_positions': [None],
            'SI': [np.nan],
            'tc_type': [np.nan]
        })
    
    # SI
    tss_list = filtered_tss['pos'].values
    total_tss = len(tss_list)
    probabilities = np.bincount(tss_list) / total_tss
    probabilities = probabilities[probabilities != 0] 
    entropy = -np.sum(probabilities * np.log2(probabilities))
    SI = 2 - entropy
    
    tc_type = 'Sharp' if SI > -1 else 'Broad'
    
    return pd.DataFrame({
        'seqnames': [seqname],
        'start': [start],
        'end': [end],
        'width': [end - start + 1],
        'total_tss': [total_tss],
        'tss_positions': [tss_list.tolist()],
        'SI': [SI],
        'tc_type': [tc_type]
    })

def main(promoters_file, tss_file, output_file):
    promoters = pd.read_csv(promoters_file, sep='\t')
    tss_positions = pd.read_csv(tss_file, sep='\t')
    
    result_list = promoters.apply(lambda row: process_promoter(row, tss_positions), axis=1).tolist()
    
    result_data = pd.concat(result_list).reset_index(drop=True)
    result_data.to_csv(output_file, sep='\t', index=False)
    
    print(result_data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process promoters and TSS positions to calculate and classify them based on Shannon Index.")
    parser.add_argument("promoters_file", help="Path to the promoters file")
    parser.add_argument("tss_file", help="Path to the TSS positions file")
    parser.add_argument("-o", "--output", required=True, help="Path to the output file")
    
    args = parser.parse_args()
    
    main(args.promoters_file, args.tss_file, args.output)
