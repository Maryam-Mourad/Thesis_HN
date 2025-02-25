nEvents=10000000
generation=0
if [ $1 ]; then
    nEvents=$1
fi

if [ $2 ]; then
    generation=$2
fi


if [ $generation = generate ]; then
    echo generating events ...
    $RAPIDSIM_ROOT/build/src/RapidSim.exe $RAPIDSIM_ROOT/HN_simulation/gamma2/gamma2 $nEvents 1
    mv gamma2_tree.root /eos/user/m/mmourad/HN_simulation/trees/
else
    echo Skipping MC generation and moving to plotting ...
fi

if [ -r plots ]; then
    echo
else
    mkdir plots
    echo plots directory created.
fi

python3 mastercode_gamma2.py $nEvents

