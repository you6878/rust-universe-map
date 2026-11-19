#!/bin/sh
# Downloads radial shells of the MEAN array of the Edenhofer et al. 2023 dust map from Zenodo.
# The array is 516 shells of 786432 big-endian floats (HEALPix nside 256, nested), from byte 5760.
# Usage: fetch_shells.sh STEP [OFFSET]   - every STEP-th shell, starting at OFFSET.
URL="https://zenodo.org/api/records/8187943/files/mean_and_std_healpix.fits/content"
SHELL_BYTES=3145728
STEP=${1:-3}; OFFSET=${2:-0}
one() {
    i=$1; f=$(printf "shells/shell_%03d.bin" $i); start=$(( 5760 + i * SHELL_BYTES ))
    tries=0
    while :; do
        have=0; [ -f $f ] && have=$(stat -f %z $f)
        [ $have -ge $SHELL_BYTES ] && return 0
        tries=$(( tries + 1 )); [ $tries -gt 100 ] && { echo "shell $i gave up"; return 1; }
        curl -sL --connect-timeout 20 -m 900 -r $(( start + have ))-$(( start + SHELL_BYTES - 1 )) "$URL" >> $f || sleep 5
    done
}
i=$OFFSET; running=0
while [ $i -lt 516 ]; do
    one $i &
    running=$(( running + 1 ))
    [ $running -ge 6 ] && { wait; running=0; }
    i=$(( i + STEP ))
done
wait
echo "DONE step $STEP offset $OFFSET: $(ls shells | wc -l | tr -d ' ') shells on disk"
