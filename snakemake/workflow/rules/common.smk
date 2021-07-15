def get_fastqs(wildcards):
    return config["samples"][wildcards.sample]


def get_salmon_outdir(wildcards, output):
    return os.path.dirname(output[0])
