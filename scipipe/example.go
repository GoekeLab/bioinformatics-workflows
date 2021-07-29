package main

import (
	"flag"
	"os"

	sp "github.com/scipipe/scipipe"
	spc "github.com/scipipe/scipipe/components"
)

func main() {
	// ----------------------------------------------------------------
	// Define input parameters to the workflow as command line flags
	// ----------------------------------------------------------------
	refPath := flag.String("ref", "", "Reference file")
	leftPath := flag.String("left", "", "The first reads file")
	rightPath := flag.String("right", "", "The second reads file")
	outDirPath := flag.String("outdir", "", "Output directory")
	flag.Parse()

	// ----------------------------------------------------------------
	// Ensure flags are set
	// ----------------------------------------------------------------
	if *refPath == "" || *leftPath == "" || *rightPath == "" || *outDirPath == "" {
		flag.PrintDefaults()
		os.Exit(1)
	}

	// ----------------------------------------------------------------
	// Initialize workflow with 4 workers
	// ----------------------------------------------------------------
	wf := sp.NewWorkflow("example-workflow", 4)

	// We need to create workflow processes to feed the files to downstream
	// processes as IPs (information packets)
	refSrc := spc.NewFileSource(wf, "ref-file", *refPath)
	leftSrc := spc.NewFileSource(wf, "left", *leftPath)
	rightSrc := spc.NewFileSource(wf, "right", *rightPath)

	// ----------------------------------------------------------------
	// Salmon indexing process
	// ----------------------------------------------------------------
	salmonIndex := wf.NewProc("salmon-index", "salmon index -t {i:ref} -i {o:index}")
	salmonIndex.In("ref").From(refSrc.Out())
	salmonIndex.SetOut("index", "data/salmon-index")

	// ----------------------------------------------------------------
	// Transcriptome alignment and quantification using Salmon
	// ----------------------------------------------------------------
	salmonAlignQuant := wf.NewProc("salmon-aligh-quant", "salmon quant -i {i:index} -l A -1 '{i:left}' -2 '{i:right}' --validateMappings -o {o:quant}")
	salmonAlignQuant.In("index").From(salmonIndex.Out("index"))
	salmonAlignQuant.In("left").From(leftSrc.Out())
	salmonAlignQuant.In("right").From(rightSrc.Out())
	salmonAlignQuant.SetOut("quant", *outDirPath)

	// ----------------------------------------------------------------
	// Quality control using FastQC
	// ----------------------------------------------------------------
	fastQC := wf.NewProc("fastqc", "mkdir -p {o:qcdir} && fastqc --quiet '{i:left}' '{i:right}' --outdir {o:qcdir}")
	fastQC.In("left").From(leftSrc.Out())
	fastQC.In("right").From(rightSrc.Out())
	fastQC.SetOut("qcdir", "data/qc")

	// ----------------------------------------------------------------
	// Run workflow
	// ----------------------------------------------------------------
	wf.Run()
}
