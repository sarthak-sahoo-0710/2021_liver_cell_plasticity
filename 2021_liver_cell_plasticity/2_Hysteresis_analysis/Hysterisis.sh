export JULIA_NUM_THREADS=11

echo "Running for reversible not Hepatocyte"
tgfb_min=25000
tgfb_max=25000
num_points=1
tgfb=$tgfb_min
temp=$((tgfb_max-tgfb_min))
dt=$((temp/num_points))
for i in $( eval echo {0..$num_points}); do 
	julia 2D_patterning_and_hysterisis_not_hep_rev.jl $tgfb out_not_hep_rev_$tgfb.txt
	echo $tgfb
	tgfb=$((tgfb + dt))
	sleep 10
done
echo "Running for Cholangiocyte"
tgfb_min=0
tgfb_max=40000
num_points=20
tgfb=$tgfb_min
temp=$((tgfb_max-tgfb_min))
dt=$((temp/num_points))
for i in $( eval echo {0..$num_points}); do 
	julia 2D_patterning_and_hysterisis_chol.jl $tgfb out_Cholangiocyte_$tgfb.txt
	echo $tgfb
	tgfb=$((tgfb + dt))
	sleep 10
done
