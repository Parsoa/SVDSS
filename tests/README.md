# (Semi-) Automatic tests

Download example data from [here](https://drive.google.com/file/d/1GD7NKz7kUhK_o2IystdCwKEBzz9r9JsB/view) and put the tar in this folder.

```
# Unzip and move input files
mkdir -p input
tar xvfz SVDSS-example.tar.gz -C input

# How to run and test a SVDSS binary
bash run-svdss.sh <path-to-SVDSS-binary> input/22.fa input/22.bam output

# Edit Dockerfile(s) - if needed (e.g., to change linux distro)

# Test binary (it must be copied here)
mv <path-to-SVDSS-binary> SVDSS_linux_x86-64
docker build -t svdss-bin:v1.0 . -f Dockerfile-bin

# Test build and then binary
docker build -t svdss-build:v1.0 . -f Dockerfile-build
```

