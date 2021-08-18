version 1.0

workflow wdl_poc {
   input {
     File reads1
     File reads2
     File ref_txome
   }
   call FastQC {
      input:
         reads1 = reads1,
         reads2 = reads2,
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

task FastQC {
  input {
     File reads1
     File reads2
  }

  command {
     mkdir fastqc_res; fastqc --quiet "${reads1}" "${reads2}" --outdir fastqc_res
  }

  output {
	 File fastqc_res = "fastqc_res"
  }

	runtime {
		docker: "biocontainers/fastqc:v0.11.9_cv8"
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

	runtime {
		docker: "salmon_docker_image_goes_here"
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

	runtime {
		docker: "salmon_docker_image_goes_here"
	}
}


