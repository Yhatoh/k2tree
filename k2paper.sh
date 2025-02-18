#!/usr/bin/env bash

echo "INDOCHINA 2004"
echo "BUILDING"
#./k2sparse.x /data/matrix/webg/indochina-2004.mtx.rowm -s 7414866
echo "MULTIPLY"
/usr/bin/time ./k2bp_mult_1.x /data/matrix/webg/indochina-2004.mtx.rowm.k2bprrr /data/matrix/webg/indochina-2004.mtx.rowm.k2bprrr
echo "--------------"

echo "ARABIC 2005"
echo "BUILDING"
#./k2bp_build.x /data/matrix/webg/arabic-2005.mtx.rowm 22744080 639999458
echo "MULTIPLY"
/usr/bin/time ./k2bp_mult_1.x /data/matrix/webg/arabic-2005.mtx.rowm.k2bprrr /data/matrix/webg/arabic-2005.mtx.rowm.k2bprrr
echo "--------------"


echo "UK 2005"
echo "BUILDING"
#./k2bp_build.x /data/matrix/webg/uk-2005.mtx.rowm 39459925 936364282
echo "MULTIPLY"
/usr/bin/time ./k2bp_mult_1.x /data/matrix/webg/uk-2005.mtx.rowm.k2bprrr /data/matrix/webg/uk-2005.mtx.rowm.k2bprrr
echo "--------------"

echo "COMPRESSED"
echo "INDOCHINA 2004"
echo "BUILDING"
#./k2sparse.x /data/matrix/webg/indochina-2004.mtx.rowm -s 7414866
echo "MULTIPLY"
/usr/bin/time ./k2bp_compr_mult.x /data/matrix/webg/indochina-2004.mtx.rowm.k2bprrri /data/matrix/webg/indochina-2004.mtx.rowm.k2bprrri
echo "--------------"

echo "ARABIC 2005"
echo "BUILDING"
./k2bp_compr.x /data/matrix/webg/arabic-2005.mtx.rowm.k2bprrr
echo "MULTIPLY"
/usr/bin/time ./k2bp_compr_mult.x /data/matrix/webg/arabic-2005.mtx.rowm.k2bprrri /data/matrix/webg/arabic-2005.mtx.rowm.k2bprrri
echo "--------------"


echo "UK 2005"
echo "BUILDING"
./k2bp_compr.x /data/matrix/webg/uk-2005.mtx.rowm.k2bprrr
echo "MULTIPLY"
/usr/bin/time ./k2bp_compr_mult.x /data/matrix/webg/uk-2005.mtx.rowm.k2bprrri /data/matrix/webg/uk-2005.mtx.rowm.k2bprrri
echo "--------------"
