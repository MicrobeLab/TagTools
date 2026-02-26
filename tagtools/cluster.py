import argparse
from sklearn.cluster import SpectralClustering
import matplotlib.pyplot as plt
import gzip
from numba import njit, prange
import random
import os
from functools import partial
import numpy as np
from collections import defaultdict, Counter
import itertools
from multiprocessing import Pool


def batched(iterable, n, *, strict=False):
    if n < 1:
        raise ValueError('n must be at least one')
    iterator = iter(iterable)
    while batch := tuple(itertools.islice(iterator, n)):
        if strict and len(batch) != n:
            raise ValueError('batched(): incomplete batch')
        yield batch


def pairwise(iterable):
    iterator = iter(iterable)
    a = next(iterator, None)
    for b in iterator:
        yield a, b
        a = b


def cmd_int(cmd):
    with os.popen(cmd) as h1:
        return int(list(h1)[0])


def uint16(data):
    return np.frombuffer(data, dtype=np.uint16)[0]


def get_header_size(bam):
    with gzip.open(bam, 'rt') as h1:
        return sum(len(i) for i in itertools.takewhile(lambda x: x.startswith('@'), h1))


def index_bgzip(file):
    total = bytearray()
    with open(file, 'rb') as h1:
        while tmp1 := h1.read(18)[-2:]:
            block_size = int(uint16(tmp1)) - 17
            tmp2 = h1.read(block_size)[-4:]
            total.extend(tmp1 + tmp2)
    total = np.frombuffer(total, dtype=[('zip', np.uint16), ('unzip', np.uint32)]).copy()
    total['zip'] += 1
    if total['zip'].sum() != os.path.getsize(file):
        raise ValueError(f'{file} file format error.')
    bam_index = np.full(len(total) * 2 + 1, len(total), dtype=np.uint64)
    bam_index[1::2] = np.cumsum(total['zip'], dtype=np.uint64)
    bam_index[2::2] = np.cumsum(total['unzip'], dtype=np.uint64)
    index_name = f'{file}.gzi'
    bam_index.tofile(index_name)


def chunk_sam(sam_file, *, chunk_num):
    start_pos = get_header_size(sam_file)
    total = int(np.fromfile(f'{sam_file}.gzi', dtype=np.uint64)[-1])
    chunk_size = int((total - start_pos) / chunk_num) + 1
    keep_pos = [start_pos]
    while start_pos < total:
        start_pos += chunk_size
        if start_pos >= total:
            keep_pos.append(total)
            break
        start_pos += cmd_int(f'bgzip -I {sam_file}.gzi -b {start_pos} {sam_file} | head -1 | wc -c')
        keep_pos.append(start_pos)
    return [(num, start, end - start) for num, (start, end) in enumerate(pairwise(keep_pos))]


def cb_acc(cb_count):
    with gzip.open(cb_count, 'rt') as h1:
        m1 = itertools.groupby((i.split() for i in h1), key=lambda x: x[0])
        return [(g_name, sum(int(i[-1]) for i in rows)) for g_name, rows in m1]


def split_by_tag(*, bam, nt, out_name, tag_name='CB', tag_list, region_file):
    mid_sam = f'{out_name}-tmp.sam.gz'
    os.system(rf'''samtools view -O BAM -@ {nt} --tag-file {tag_name}:{tag_list} -L {region_file} {bam} | \
        samtools sort -@ {nt} -t {tag_name} -O SAM -T {mid_sam}-tmp -o {mid_sam}''')
    index_bgzip(mid_sam)
    regions = chunk_sam(sam_file=mid_sam, chunk_num=nt * 2)
    header = f'{mid_sam}.header'
    os.system(rf'''samtools view -H {bam} > {header}''')
    
    cmd1 = [rf'''bgzip -b {start} -s {q_len} {mid_sam} | \
        mawk -F '\t' '{{
            for(i=12;i<=NF;i++){{
                if(substr($i, 1, 2)=="{tag_name}"){{ tag=substr($i, 6); print tag"\t"length+1; next;}}
            }} }}' | bgzip > {mid_sam}-{num}.wc.gz ''' for num, start, q_len in regions]
    
    header_size = get_header_size(mid_sam)
    with Pool(nt) as p1:
        p1.map(os.system, cmd1, chunksize=1)
        cb_result_files = [f'{mid_sam}-{num}.wc.gz' for num, *_ in regions]
        final = list(itertools.chain.from_iterable(p1.map(cb_acc, cb_result_files, chunksize=1)))
        cell_order = []
        cell_count = defaultdict(int)
        for cell, count in final:
            cell_count[cell] += count
            if cell not in cell_order:
                cell_order.append(cell)
        cell_len = [cell_count[i] for i in cell_order]
        c_index = [(start, end - start) for start, end in pairwise(itertools.accumulate(cell_len, initial=header_size))]
        c_index = [(cell, start, q_len) for cell, (start, q_len) in zip(cell_order, c_index)]
        cmd2 = [rf'''bgzip -b {start} -s {q_len} {mid_sam} | \
            mawk -F '\t' 'BEGIN{{OFS="\t"}} {{ if($5==255){{$5=60;}} print;}}' | cat {header} - | \
            samtools view -O BAM -o {mid_sam}-{cell}.bam --write-index ''' for cell, start, q_len in c_index]
        p1.map(os.system, cmd2, chunksize=1)
    os.system(rf'''rm {' '.join(cb_result_files)} {mid_sam} {mid_sam}.gzi {header}''')
    return [(cell, f'{mid_sam}-{cell}.bam') for cell in cell_order]


def call_pl(*, bam, nt, out_name, ref, tag_name, tag_list, site):
    tmp_dir = f'{out_name}-dir'
    if os.path.exists(tmp_dir):
        raise FileExistsError(f'{tmp_dir} exits.')
    os.system(f'mkdir {tmp_dir}')
    split_name = f'{tmp_dir}/for-split'
    
    new_size = f'{split_name}.site.gz'
    os.system(rf'''samtools depth -@ {nt} -b {site} {bam} | mawk '$NF>0{{print $1"\t"$2}}' | \
            tabix -R - {site} | mawk '!a[$1" "$2]++' | bgzip > {new_size} && tabix -s1 -b2 -e2 -C {new_size}
            gzip -dc {new_size} | sort -k1,1 -k2,2n | mawk '{{print NR"\t"$0}}' | bgzip > {new_size}-num.gz
            tabix -s1 -b2 -e2 -C {new_size}-num.gz''')
    
    cell_bam = split_by_tag(bam=bam, nt=nt, out_name=split_name, tag_name=tag_name, tag_list=tag_list,
                            region_file=new_size)
    random.shuffle(cell_bam)
    b_size = 20
    keep_cell = list(batched([i[0] for i in cell_bam], b_size))
    keep_bam = list(batched([i[1] for i in cell_bam], b_size))
    
    call_pl_cmds = [rf'''samtools depth -b {new_size} {' '.join(bam_list)} | \
        mawk '{{for (i=3; i<=NF; i++) if($i > 0){{print $1"\t"$2; break}} }}' | tabix -R - {new_size} | \
        bgzip > {split_name}-chunk{num}.site.gz && tabix -s1 -b2 -e2 -C {split_name}-chunk{num}.site.gz
        bcftools mpileup --ignore-RG -f {ref} -I -E -T {split_name}-chunk{num}.site.gz -Ob {' '.join(bam_list)} | \
        bcftools call -Am -C alleles -T {split_name}-chunk{num}.site.gz -Ob | \
        bcftools view -e 'INFO/DP<1' -Ob | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t[%GT:%PL\t]\n" - | \
        mawk '
            BEGIN {{
                sample_ids_str = "{' '.join(cell_list)}"
                split(sample_ids_str, sample_ids, " ")
                sample_count = length(sample_ids)
            }}
            {{
                for (i = 1; i <= sample_count; i++) {{
                    pl = $(i+4)
                    if(substr(pl,1,1) !="."){{
                        split(substr(pl,5), nums, ",")
                        print sample_ids[i]"\t"$1"\t"$2"\t"$3"\t"$4"\t"nums[1]"\t"nums[2]"\t"nums[3]
                    }}
                }}
            }}' | sort -k1,1 -k2,2 -k3,3n | mawk '{{print NR"\t"$0}}' | bgzip > {split_name}-chunk{num}.col.pl.gz
        rm {split_name}-chunk{num}.site.gz {split_name}-chunk{num}.site.gz.csi
        rm {' '.join(bam_list)} {' '.join(f'{i}.csi' for i in bam_list)}
        ''' for num, (cell_list, bam_list) in enumerate(zip(keep_cell, keep_bam))]
    
    keep_pl = [rf'''{split_name}-chunk{num}.col.pl.gz''' for num, _ in enumerate(keep_cell)]
    with Pool(nt) as p1:
        p1.map(os.system, call_pl_cmds, chunksize=1)
    pl_file = f'{out_name}.pl.gz'
    os.system(f'rm {new_size} {new_size}.csi')
    os.system(rf'''cat {' '.join(keep_pl)}  > {pl_file} && tabix -s2 -b1 -e1 -C {pl_file}''')
    os.system(f'rm -R {tmp_dir}')
    return pl_file


def chunk_to_numpy(cell, *, pos_file, pl_file):
    with os.popen(f'tabix {pl_file} {cell}') as h1:
        gl = 10 ** (np.array([i.split()[-3:] for i in h1], dtype=np.uint32) * -0.1)
        gp = gl / (gl.sum(axis=1)[:, np.newaxis])  # # convert GL to GP
        ds = (gp * np.array([0, 1, 2])).sum(axis=1)  # covert GP to DS
    
    cmd1 = rf'''tabix {pl_file} {cell}|mawk '{{print $3"\t"$4}}'|tabix -R - {pos_file} |perl -nale 'print $F[-1]' '''
    with os.popen(cmd1) as h1:
        pos_nums = np.array([i for i in h1], dtype=np.uint32)
    
    d1 = np.full(len(ds), 0, dtype=np.dtype([('pos', np.uint32), ('ds', np.float64)]))
    d1['pos'] = pos_nums
    d1['ds'] = ds
    return d1


def pl_to_numpy(*, pl_file, nt, pos_file):
    with os.popen(f'tabix -l {pl_file}') as h1:
        all_cell = [i.strip() for i in h1]
    fun1 = partial(chunk_to_numpy, pl_file=pl_file, pos_file=pos_file)
    keep_min_snp = 50
    with Pool(nt) as p1:
        final = p1.map(fun1, all_cell, chunksize=10)
        keep_len = [len(i) for i in final]
        m1 = [i > keep_min_snp for i in keep_len]
        cell_order = list(itertools.compress(all_cell, m1))
        final = list(itertools.compress(final, m1))
        keep_len = list(itertools.compress(keep_len, m1))
    
    pl_bin = np.concatenate(final)
    acc = itertools.accumulate(keep_len, initial=0)
    region_data = np.fromiter(itertools.pairwise(acc), dtype=[(it1, np.uint32) for it1 in 'start end'.split()])
    
    index_data = get_combin(len(region_data))
    a = index_data[::2]
    b = index_data[1::2]
    index_data = np.zeros(len(a), dtype=[(it1, np.uint32) for it1 in 'index1 index2'.split()])
    index_data['index1'] = a
    index_data['index2'] = b
    
    cor_matrix = get_half_cor_matrix(pl_bin=pl_bin, region_data=region_data, index_data=index_data)
    cor_matrix = get_full_cor_matrix(cor_matrix)
    return cor_matrix, cell_order


@njit
def cosine_similarity(a, b):
    dot_product = 0.0
    norm_a = 0.0
    norm_b = 0.0
    for i in range(len(a)):
        dot_product += a[i] * b[i]
        norm_a += a[i] ** 2
        norm_b += b[i] ** 2
    if norm_a == 0.0 or norm_b == 0.0:
        return 0.0
    else:
        return dot_product / (np.sqrt(norm_a) * np.sqrt(norm_b))


def get_full_cor_matrix(final):
    cell_num = int(np.sqrt(2 * len(final)))
    full_matrix = np.zeros((cell_num, cell_num), dtype=np.float64)  # empty N X N matrix
    full_matrix[np.triu_indices(cell_num)] = final  # fill the up half matrix
    full_matrix += full_matrix.T - np.diag(full_matrix.diagonal())  # fill the bot half matrix
    return full_matrix


@njit(parallel=True)
def get_half_cor_matrix(*, pl_bin, region_data, index_data):
    final = np.zeros(len(index_data), dtype=np.float64)
    for num1 in prange(len(index_data)):
        index1 = index_data[num1]['index1']
        index2 = index_data[num1]['index2']
        s1 = region_data[index1]['start']
        e1 = region_data[index1]['end']
        s2 = region_data[index2]['start']
        e2 = region_data[index2]['end']
        tmp1 = pl_bin[s1:e1]
        tmp2 = pl_bin[s2:e2]
        share1, share2 = get_share(tmp1, tmp2)
        final[num1] = cosine_similarity(share1['ds'], share2['ds'])
    return final


@njit
def get_combin(num_size):
    cur_pos = num_size
    total = 0
    while cur_pos > 0:
        total += cur_pos
        cur_pos -= 1
    tmp1 = np.zeros(total * 2, dtype=np.uint32)
    start = 0
    index = 0
    while start < num_size:
        for num2 in range(start, num_size):
            tmp1[index] = start
            index += 1
            tmp1[index] = num2
            index += 1
        start += 1
    return tmp1


@njit()
def get_share(num1, num2):
    tmp1 = np.zeros(len(num1) + len(num2), num1.dtype)
    tmp2 = np.zeros(len(num1) + len(num2), num1.dtype)
    index1 = 0
    index2 = 0
    keep_count = 0
    
    while index1 < len(num1) and index2 < len(num2):
        row1 = num1[index1]
        row2 = num2[index2]
        
        if row1['pos'] == row2['pos']:
            tmp1[keep_count] = row1
            tmp2[keep_count] = row2
            index1 += 1
            index2 += 1
            keep_count += 1
        elif row1['pos'] < row2['pos']:
            index1 += 1
        else:
            index2 += 1
    return tmp1[:keep_count], tmp2[:keep_count]


@njit(parallel=True)
def predict(*, big_matrix, small_data):
    want_index = np.zeros(len(big_matrix), dtype=np.uint32)
    for num in prange(len(big_matrix)):
        cell_val = big_matrix[num]  # select the col to predict the nearest point
        tmp_row = np.zeros(len(small_data), dtype=np.float64)
        for num2 in range(len(small_data)):
            distance = cosine_similarity(cell_val, small_data[num2])
            tmp_row[num2] = distance
        keep_index = sorted([(cor_val, index) for index, cor_val in enumerate(tmp_row)], reverse=True)[0][-1]
        want_index[num] = keep_index
    return want_index


def cluster_predict(*, cor_matrix, cluster_num):
    # replace the clustering method with SpectralClustering
    spectral = SpectralClustering(n_clusters=cluster_num, affinity='precomputed', random_state=42)
    labels = spectral.fit_predict(cor_matrix)
    centers = np.array([cor_matrix[labels == i].mean(axis=0) for i in range(cluster_num)])
    combine_cell = get_combin(cluster_num)
    combine_cell.shape = (-1, 2)
    new_data = np.array([(centers[a] + centers[b]) / 2 for a, b in combine_cell])
    keep_index = predict(big_matrix=cor_matrix, small_data=new_data)
    final = combine_cell[keep_index]
    return final  # return the cell combine index tuple (1,2) (1,1) ...


def make_new_centers(*, data, final):
    m1 = [(a, b) for a, b in final]
    g_names = sorted({(a, b) for a, b in m1 if a == b})
    groups = []  # mask group-columns, select the non-double cells
    for i in g_names:
        m2 = np.array([i == j for j in m1], dtype=np.bool_)
        groups.append(m2)
    
    new_center = []
    for i in groups:
        sub1 = data[i]
        new_center.append(sub1.mean(axis=0))
    combine_cell = get_combin(len(new_center))
    combine_cell.shape = (-1, 2)
    new_center = np.concatenate([(new_center[a] + new_center[b]) / 2 for a, b in combine_cell])
    new_center.shape = (len(combine_cell), -1)
    keep_index = predict(big_matrix=data, small_data=new_center)
    new_final = combine_cell[keep_index]
    return new_final


def get_predict_index(*, out_name, cluster_num, pl_file, nt, plot_log=False):
    pos_file = f'{out_name}.pos.gz'
    os.system(rf'''bgzip -dc {pl_file} | mawk '!a[$3" "$4]++{{print $3"\t"$4}}'| \
            sort -k1,1 -k2,2n | perl -nale 'print "$_\t$."' | bgzip > {pos_file}
            tabix -s1 -b2 -e2 -C {pos_file}''')
    
    matrix_data, cell_order = pl_to_numpy(pl_file=pl_file, nt=nt, pos_file=pos_file)
    tmp_final = cluster_predict(cor_matrix=matrix_data, cluster_num=cluster_num)
    new_index = make_new_centers(data=matrix_data, final=tmp_final)
    new_index = make_new_centers(data=matrix_data, final=new_index)
    with open(f'{out_name}-predict.txt', 'wt') as h1:
        for a, b in zip(cell_order, new_index):
            print(a, *b, sep='\t', file=h1)
    os.system(rf'''rm {pos_file} {pos_file}.csi
    perl -nale 'print join "\t",($F[0], "group$F[-1]") if $F[-2] eq $F[-1]' \
    {out_name}-predict.txt > {out_name}-predict.singlet''')
    
    if plot_log:
        # Elbow Method plot
        inertia = []
        k_values = np.array(list(range(2, 11)), dtype=np.uint32)
        for k in k_values:
            spectral = SpectralClustering(n_clusters=k, affinity='precomputed', random_state=42)
            labels = spectral.fit_predict(matrix_data)
            cluster_centers = np.array([matrix_data[labels == i].mean(axis=0) for i in range(k)])
            inertia_mtp = 0
            for i in range(k):
                cluster_data = matrix_data[labels == i]
                distances = np.linalg.norm(cluster_data - cluster_centers[i], axis=1)
                inertia_mtp += np.sum(distances ** 2)
            inertia.append(inertia_mtp)
        
        inertia = np.array(inertia, dtype=np.float64)
        inertia = (inertia - np.min(inertia)) / (np.max(inertia) - np.min(inertia))  # Min-Max normalized
        singlet_num = [a for a, b in new_index if a == b]
        cluster_count = Counter(singlet_num)
        names = list(cluster_count.keys())
        counts = list(cluster_count.values())
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
        ax1.plot(k_values, inertia, marker='o')
        ax1.set_xlabel('Number of clusters (k)')
        ax1.set_ylabel('Inertia')
        ax1.set_title('Elbow Method For Optimal k')
        ax2.pie(counts, labels=names, autopct='%1.1f%%', startangle=140)
        ax2.set_title(f'Distribution of cluster (n={len(singlet_num)})')
        plt.tight_layout()
        plt.savefig(f'{out_name}.inertia.pdf')
        plt.close()
        
        with open(f'{out_name}-best-k.txt', 'wt') as h1:
            for num1, num2 in zip(k_values, inertia):
                print(out_name, num1, num2, file=h1)
    return out_name


def run_cluster_predict(*, bam, out_name, cb_list, site, nt, ref, sample_num, plot_log=False):
    pl_file = call_pl(bam=bam, out_name=out_name, site=site, ref=ref, nt=48, tag_name='CB', tag_list=cb_list)
    sample_num = int(sample_num)
    get_predict_index(pl_file=pl_file, out_name=out_name, cluster_num=sample_num, nt=nt, plot_log=plot_log)


def main(args=None):
    """Main entry point. If args is None, parse from command line."""
    parser = argparse.ArgumentParser(description='Cluster prediction from SNP data')
    # 必需参数（没有默认值，required=True）
    parser.add_argument('--bam', required=True, help='BAM file, should be sorted and index by samtools')
    parser.add_argument('--out-name', required=True, help='out file name')
    parser.add_argument('--cb-list', required=True, help='cell barcode list file')
    parser.add_argument('--snp-site', required=True, help='common SNPs site file, indexed with tabix')
    parser.add_argument('--ref', required=True, help='reference genome file')
    parser.add_argument('--sample-num', type=int, required=True, help='number of sample multiplexed')
    # 可选参数（有默认值，required=False）
    parser.add_argument('--nt', '--threads', type=int, default=16, help='number of threads to used, default: 16')
    parser.add_argument('--plot-log', action='store_true', help='plot the log file, default: false')
    
    if args is None:
        args = parser.parse_args()
    # 调用主要函数
    return run_cluster_predict(
        bam=args.bam,
        nt=args.nt,
        out_name=args.out_name,
        site=args.snp_site,
        ref=args.ref,
        sample_num=args.sample_num,
        plot_log=args.plot_log,
        cb_list=args.cb_list
    )


if __name__ == '__main__':
    main()
