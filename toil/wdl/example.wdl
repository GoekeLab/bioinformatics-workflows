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
}


task SalmonIndex {
  input {
     File ref_txome
  }

  command {
     salmon index -t "${ref_txome}" -i index
  }

  output {
	 File index = "index"
  }
}

task SalmonAlignQuant {
  input {
	 File reads1
	 File reads2
     File index
  }

  command {
     salmon quant -i "${index}" -l A -1 "${reads1}" -2 "${reads2}" --validateMappings -o quant
  }

  output {
	 File quant = "quant"
  }
}


