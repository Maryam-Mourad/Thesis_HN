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
    $RAPIDSIM_ROOT/build/src/RapidSim.exe $RAPIDSIM_ROOT/HN_simulation/gamma3/gamma3 $nEvents 1
    mv gamma3_tree.root /eos/user/m/mmourad/HN_simulation/trees/
else
    echo Skipping MC generation and moving to plotting ...
fi
python3 mastercode_gamma3.py $nEvents

