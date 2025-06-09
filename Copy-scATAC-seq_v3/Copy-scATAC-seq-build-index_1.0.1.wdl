version 1.0

workflow main{

	input{
		String refName
		File GTF
		String chrM="chrM"
		File FASTA
		File blacklist="None"
		Int Cpu = 4
		Int Mem = 100
	}
  String dockerUrl="public-library/tianyongjing_bbced6197e9a469f9b9d4a3c9006a400_public:latest"

  call checkref{
		input:
		gtf = GTF,
    fasta=FASTA,
		dockerUrl = dockerUrl,
		cpu = Cpu,
		mem = Mem
	}

	call buildIndex{
		input:
		refcode = refName,
		gtf = GTF,
		chrm = chrM,
		fasta = FASTA,
		blacklist = blacklist,
    refstatistic=checkref.refstatistic,
    dockerUrl=dockerUrl,
		cpu = Cpu,
		mem = Mem
	}
	output{
		File refreport=checkref.refreport
		File refstatistic=checkref.refstatistic
		File star_index = buildIndex.index_dir	
	}
}

task checkref{
    input {
		String gtf
    String fasta
		String outdir="./"
		String dockerUrl
		Int cpu
		Int mem
}
command<<<
    gtf2gff ~{gtf} > ~{outdir}/checkgtf.gff
    gff3_QC --gff ~{outdir}/checkgtf.gff --fasta ~{fasta} --output ~{outdir}/report.txt --statistic ~{outdir}/statistic.txt
    awk '
    {
      if ($3 == "Error" && $1 != "Esf0014") {
        print "Error found in line " NR ": " $0
        exit 1
      }
    }
    ' ~{outdir}/statistic.txt 1>&2
>>>
	runtime{
		docker_url:"${dockerUrl}"
		req_cpu:cpu
		req_memory:"${mem}Gi"
	}
		output{
		File refreport="${outdir}/report.txt"
		File refstatistic="${outdir}/statistic.txt"
	}
}

task buildIndex{
		input {
		String refcode
		String gtf
		String fasta
		String chrm
		String blacklist
    String dockerUrl
    File refstatistic
    String outdir="./"
		Int cpu
		Int mem
		}
		command<<<
		mkdir -p ~{outdir}/scATAC_ref
		if [ "~{blacklist}" = "None" ]
		then
		echo "None"
		dnbc4tools atac mkref --ingtf ~{gtf} --fasta ~{fasta}  --genomeDir ~{outdir}/scATAC_ref --species "~{refcode}" --chrM ~{chrm} && \
		cp ~{fasta} ~{outdir}/scATAC_ref/genome.fasta && \
		cp ~{gtf}  ~{outdir}/scATAC_ref/genes.gtf		
		else
		echo "Blacklist"
		dnbc4tools atac mkref --ingtf ~{gtf} --fasta ~{fasta}  --genomeDir ~{outdir}/scATAC_ref --species "~{refcode}" --chrM ~{chrm} --blacklist ~{blacklist} && \
		cp ~{fasta} ~{outdir}/scATAC_ref/genome.fasta && \
		cp ~{gtf}  ~{outdir}/scATAC_ref/genes.gtf && \
		cp ~{blacklist}  ~{outdir}/scATAC_ref/blacklist.bed
		fi
		>>>
	runtime{
		docker_url: "${dockerUrl}"
		req_cpu:cpu
		req_memory:"${mem}Gi"
	}
		output{
		File index_dir="${outdir}/scATAC_ref"
	}
}