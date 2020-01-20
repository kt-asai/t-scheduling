#!/bin/bash
if [ ! -d opt-tpar-gauss-false-2 ]; then
        mkdir opt-tpar-gauss-false-2
fi
for f in `find ./*.qc`
do
  ../build/t-scheduling $f tpar gauss false 2 > opt-tpar-gauss-false-2/$f.opt
done

if [ ! -d opt-tpar-gauss-false-4 ]; then
        mkdir opt-tpar-gauss-false-4
fi
for f in `find ./*.qc`
do
  ../build/t-scheduling $f tpar gauss false 4 > opt-tpar-gauss-false-4/$f.opt
done

if [ ! -d opt-tpar-gauss-false-6 ]; then
        mkdir opt-tpar-gauss-false-6
fi
for f in `find ./*.qc`
do
  ../build/t-scheduling $f tpar gauss false 6 > opt-tpar-gauss-false-6/$f.opt
done

if [ ! -d opt-tpar-par-true-2 ]; then
        mkdir opt-tpar-par-true-2
fi
for f in `find ./*.qc`
do
  ../build/t-scheduling $f tpar par true 2 > opt-tpar-par-true-2/$f.opt
done

if [ ! -d opt-tpar-par-true-4 ]; then
        mkdir opt-tpar-par-true-4
fi
for f in `find ./*.qc`
do
  ../build/t-scheduling $f tpar par true 4 > opt-tpar-par-true-4/$f.opt
done

if [ ! -d opt-tpar-par-true-6 ]; then
        mkdir opt-tpar-par-true-6
fi
for f in `find ./*.qc`
do
  ../build/t-scheduling $f tpar par true 6 > opt-tpar-par-true-6/$f.opt
done

# if [ ! -d opt-tskd-par-true-2 ]; then
#         mkdir opt-tskd-par-true-2
# fi
# for f in `find ./*.qc`
# do
#   ../build/t-scheduling $f tskd par true 2 > opt-tskd-par-true-2/$f.opt
# done

# if [ ! -d opt-tskd-par-true-4 ]; then
#         mkdir opt-tskd-par-true-4
# fi
# for f in `find ./*.qc`
# do
#   ../build/t-scheduling $f tskd par true 4 > opt-tskd-par-true-4/$f.opt
# done

# if [ ! -d opt-tskd-par-true-6 ]; then
#         mkdir opt-tskd-par-true-6
# fi
# for f in `find ./*.qc`
# do
#   ../build/t-scheduling $f tskd par true 6 > opt-tskd-par-true-6/$f.opt
# done