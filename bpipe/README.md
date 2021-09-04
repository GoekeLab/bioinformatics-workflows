# Example Bpipe Workflow

[Bpipe](http://bpipe.org) is a workflow manager focused on emulating the simplicity of 
the command line and making your pipelines look as close to the literal commands 
that run as possible.

Bpipe pipeline are written in Groovy, a scripting version of Java that enables powerful
yet high performing DSLs to be easily developed. The definition of the workflow is 
in the [quantify_rna.groovy](quantify_rna.groovy) file.

To run this workflow, you can first install Bpipe using [SDKMan](https://sdkman.io/sdks#bpipe) 
and then execute from within the bpipe directory:

```
bpipe run quantify_rna.groovy  ../test_data/*.fq.gz
```

Please note the default configuration will run the commands on the local computer. It is easy
to configure Bpipe for other environments, and an example is shown in the 
[bpipe.config](bpipe.config) file which illustrates how to configure it to run using
PBS Torque. To run the Torque version you can tell Bpipe to use the alternative environment like
so:

```
bpipe run --env torque  quantify_rna.groovy  ../test_data/*.fq.gz
```

If you would like to see how the HTML report that Bpipe makes looks, add the `-r` option:

```
bpipe run -r --env torque  quantify_rna.groovy  ../test_data/*.fq.gz
```

This example only uses very simple features of Bpipe. Bpipe has many more features
which you can explore in the [documentation](http://docs.bpipe/org).

Other useful commands to experiment with are:

- `bpipe log` (see logs of a running pipeline)
- `bpipe stop` (stop a running pipeline)
- `bpipe history` (show history)


Thanks for trying out the Bpipe example!


