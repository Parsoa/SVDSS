# Portable binary

To build a portable binary that works on any linux machine (at least we hope so),
we use the [holy-build-box](https://github.com/phusion/holy-build-box) system.

```
docker pull phusion/holy-build-box-64:latest
# cd root of this repo
docker run -t -i --rm -v `pwd`:/io phusion/holy-build-box-64:3.0.4 /hbb_exe/activate-exec bash /io/holy/holy-compile.sh
# quite useful (i.e., to run libcheck on a binary):
# docker run -t -i --rm -v `pwd`:/io phusion/holy-build-box-64:3.0.4 /hbb_exe/activate-exec bash
```

The binary is then copied here (`holy` subfolder).

Note:
* this binary should not be linked to any [non-essential library](https://github.com/phusion/holy-build-box/blob/master/ESSENTIAL-SYSTEM-LIBRARIES.md) (at least this is what `libcheck` says)
* to do so, we had to "manually" compile libbz2 and liblzma (see the workaround in [holy-compile.sh](holy-compile.sh))
* we then modified `CMakeLists.txt` accordingly (see `HOLYBUILD` option)