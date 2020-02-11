import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-k", help="kaiju output file")
parser.add_argument("-s", help="query fastq/fasta file used by kaiju")
parser.add_argument("-o", help="path name for output fasta file of extracted sequences")
args = parser.parse_args()

mapped_reads = {i.strip().split('\t')[1]:0 for i in open(args.k)}

print "Finished processing kaiju output file!!!"

read_count = 0
av_read_len = 0

with open(args.o,'w') as out, open('read_file_info.txt','w') as out2:
	if args.s.endswith('fastq'):
		for h,i in enumerate(open(args.s)):
			if h%4 == 0:
				id = i[1:].split(' ')[0].strip()
			elif h%4 == 1:
				s = i.strip()
				av_read_len += len(s)
				read_count += 1.0
				if id in mapped_reads:
					out.write(">"+id+"\n"+s+"\n")
		out2.write(str(read_count)+'\n'+str(av_read_len/read_count))
	elif args.s.endswith('fasta'):
		from Bio import SeqIO
		for i in SeqIO.parse(args.s,'fasta'):
			id = str(i.id)
			s = str(i.seq)
			read_count += 1
			av_read_len += len(s)
			if id in mapped_reads:
				out.write(">"+id+"\n"+s+"\n")
		out2.write(str(read_count)+'\n'+str(av_read_len/read_count))
