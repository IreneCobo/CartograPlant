#Screen command
screen #opens a new screen session
screen -S name-of-session #put name to the new session
#Ctrl+a+d to leave a session
screen -r/screen -ls #check preexisting sessions
screen -r name-of-session #come back to a preexisting session

#Download SRA sequences from NCBI in FASTQ (Panel1-3, accesion number SRP055524, 411 sequences, SRR1819206-SRR1819616)
module avail
module load sratoolkit/2.8.2
module display sratoolkit/2.8.2 #to get the program PATH, in my case /isg/shared/apps/sratoolkit/2.8.2/bin/
for ((i=206;i<=616;i++))
do 
/isg/shared/apps/sratoolkit/2.8.2/bin/fastq-dump --split-files SRR1819$i
done

#Compress fastq sequences to fastq.gz
gzip *fastq

#


