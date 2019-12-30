import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-k", help="kaiju output file")
parser.add_argument("-s", help="query fastq/fasta file used by kaiju")
parser.add_argument("-o", help="path name for output fasta file of extracted sequences")
args = parser.parse_args()

mapped_reads = {i.strip().split('\t')[1]:0 for i in open(args.k)}

print "Finished processing kaiju output file!!!"

with open(args.o,'w') as out:
	if args.s.endswith('fastq'):
		for h,i in enumerate(open(args.s)):
			if h%4 == 0:
				id = i[1:].split(' ')[0].strip()
			elif h%4 == 1:
				s = i.strip()
				if id in mapped_reads:
					out.write(">"+id+"\n"+s+"\n")
	elif args.s.endswith('fasta'):
		from Bio import SeqIO
		for i in SeqIO.parse(args.s,'fasta'):
			id = str(i.id)
			if id in mapped_reads:
				out.write(">"+str(i.id)+"\n"+str(i.seq)+"\n")
