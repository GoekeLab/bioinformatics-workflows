version 1.0

workflow wdl_poc {
   input {
     File reads1
     File reads2
     File ref_txome
   }
   call FastQCone {
      input:
         reads = reads1,
   }
   call FastQCtwo {
      input:
         reads = reads2,
   }
   call SalmonIndex {
     input:
        ref_txome = ref_txome,
  }
   call SalmonAlignQuant {
     input:
        reads1 = reads1,
		reads2 = reads2,
		index = SalmonIndex.index
  }
}

task FastQCone {
  input {
     File reads
  }

  command {
     zcat "${reads}" | fastqc stdin:readsone
  }

  output {
	 File fastqc_res = "readsone_fastqc.html"
  }
  
  runtime {
     docker: 'pegi3s/fastqc'
  }
}

task FastQCtwo {
  input {
     File reads
  }

  command {
     zcat "${reads}" | fastqc stdin:readstwo
  }

  output {
	 File fastqc_res = "readstwo_fastqc.html"
  }

  runtime {
     docker: 'pegi3s/fastqc'
  }
}


task SalmonIndex {
  input {
     File ref_txome
     String index = "index"
  }

  command {
     salmon index -t "${ref_txome}" -i ${index}; tar -cvzf index.tar.gz ${index}
  }

  output {
	 File index = "index.tar.gz"
  }
  
  runtime {
     docker: 'combinelab/salmon'
  }
}

task SalmonAlignQuant {
  input {
	 File reads1
	 File reads2
     File index
     String indexdir = "index"
     String quantex = "quant"
  }

  command {
     tar -xzf ${index};
     salmon quant -i "${indexdir}" -l A -1 "${reads1}" -2 "${reads2}" --validateMappings -o quant;
     tar -cvzf quant.tar.gz ${quantex}
  }

  output {
	 File quant = "quant.tar.gz"
  }

  runtime {
     docker: 'combinelab/salmon'
  }
}


