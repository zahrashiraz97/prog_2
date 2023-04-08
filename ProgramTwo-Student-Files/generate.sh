#! /bin/bash

# a "wrapper script" for the nextgen program that accomplishes two things
# - fixes the number of edges in the header for geometric graphs
# - creates an output file with a standardized name
#   of the form tt_nnnnn_mmmmm_ss_wwww.gph
#   where tt is the type (rn, gm, or gw)
#         nnnnn is the number of nodes (with leading 0's if less than 10000)
#         mmmmm is the number of edges
#         ss is the seed
#         wwww is the maximum edge weight (or 0000 if not specified)
# If really large graphs are in use, the script can be edited to modify the number of digits
# in relevant parts of the name.
# ASSUME that the nextgen program is in the same directory as this script
#                 and normalize_weights as well if MAX_WEIGHT argument is used

# max number of attempts is MAX_ATTEMPT_MULTIPLIER * number of instances requested
MAX_ATTEMPT_MULTIPLIER=3
# minimum number of attempts to achieve number of instances requested; overrides the multiplier
MIN_ATTEMPTS=50
# maximum 32 or 64 bit integer, respectively
MAX_32_BIT=2147483647
MAX_64_BIT=9223372036854775807

usage() {
    echo "$1 [ NUMBER_OF_INSTANCES ] TYPE NODES EDGES STARTING_SEED [ MAX_WEIGHT ]"
    echo " creates random graphs with the given number of nodes and edges"
    echo " NUMBER_OF_INSTANCES is the number of graphs created, an integer (=1 if absent)"
    echo " TYPE is one of:"
    echo "  rn (random), gm (geometric), gw (geometric with wraparound), or dt (Delaunay triangulation)"
    echo " NODES and EDGES are number of nodes and edges, respectively"
    echo " STARTING_SEED is the random seed of the first (attempted) instance"
    echo "               remaining seeds are consecutive integers after that"
    echo " MAX_WEIGHT (optional) is the desired maximum edge weight"
    echo "            -1 or -2 means that the sum of n-1 edges (n = NODES)"
    echo "            will not exceed the maximum 32 or 64 bit integer, respectively"
    echo "*** Individual graph generation may fail if a geometric graph is not connected;"
    echo "    in that case, there is an appropriate error message."
    echo "    The script keeps trying until enough graphs are generated or it decides to give up."
    echo "    Current limits are $MAX_ATTEMPT_MULTIPLIER * NUMBER_OF_INSTANCES"
    echo "    or $MIN_ATTEMPTS, whichever comes *last*."
}

if [ $# -lt 4 ] || [ $# -gt 6 ]; then
    usage $0
    exit 1
fi

if [[ $1 =~ [0-9] ]]; then
    num_instances=$1
    shift
else
    num_instances=1
fi

graph_type=$1
shift
num_nodes=$1
node_string=`echo $num_nodes | awk '{printf "%05d",$1}'`
shift
num_edges=$1
shift
starting_seed=$1
if [ $starting_seed -lt 0 ]; then starting_seed=$((-starting_seed)); fi
shift
if [ $# -ne 0 ]; then
    max_weight=$1
    weight_string=`echo $max_weight | awk '{printf "%04d",$1}'`
    unnormalized_file=yy-$$.gph
    if [ $max_weight = '-1' ]; then
        max_weight=`echo "$MAX_32_BIT / ($num_nodes-1)" | bc`
        weight_string="32bt"
    elif [ $max_weight = '-2' ]; then
        max_weight=`echo "$MAX_64_BIT / ($num_nodes-1)" | bc`
        weight_string="64bt"
    elif [ $max_weight -lt 0 ]; then
        echo "bad maximum weight $max_weight"
        usage "./generate.sh"
        exit 1
    fi
else
    max_weight=0
    weight_string=`echo $max_weight | awk '{printf "%04d",$1}'`
fi

temp_file=zz-$$.gph

# make sure the nextgen is gotten from the script directory,
# in case the script is called from elsewhere
script_directory=${0%/*}
executable=$script_directory/nextgen
normalize=$script_directory/normalize-weights

# convert abbreviation for graph type into one that nextgen understands
case $graph_type in
    rn | random)
        gen_type=random
        ;;
    gm | geometric)
        gen_type=geometric
        ;;
    gw | geo_wrap)
        gen_type=geo_wrap
        ;;
    dt | delaunay | triangulation)
        echo "Warning: number of edges ignored in Delaunay triangulation"
        executable=$script_directory/triangulation
        gen_type=""
        num_edges=""
        ;;
    *)
        echo "unknown graph type $graph_type"
        usage $0
        exit 1
        ;;
esac

max_attempts=$((MAX_ATTEMPT_MULTIPLIER * $num_instances))
if [ $max_attempts -lt $MIN_ATTEMPTS ]; then
    max_attempts=$MIN_ATTEMPTS
fi

attempts=0
successful_attempts=0
seed=$starting_seed
until [[ $successful_attempts -ge $num_instances ]] || [[ $attempts -gt $max_attempts ]]; do
    attempts=$((attempts + 1))
    edge_string=`echo $num_edges | awk '{printf "%06d",$1}'`
    seed_string=`echo $seed | awk '{printf "%02d",$1}'`
    $executable $gen_type $num_nodes $num_edges $seed > $temp_file
    success=$?
    if [ $success -gt 1 ]; then
        # nextgen returned a value that indicates no attempt is going to succeed
        echo "*** Aborted - see message from program ***"
        exit $success
    fi
    if [ $success -eq 0 ] && [ $max_weight -gt 0 ]; then
        mv $temp_file $unnormalized_file
        $normalize $max_weight < $unnormalized_file > $temp_file
        rm $unnormalized_file
    fi
    if [ $success -eq 0 ]; then
        out_file=${graph_type}_${node_string}_${edge_string}_${weight_string}_${seed_string}.gph
        if [ $graph_type = "random" ]; then
            mv $temp_file $out_file
        else
            # geometric or triangulation: need to fix number of edges and revise file name
            # count lines that begin with e followed by a space
            desired_num_edges=$num_edges
            true_num_edges=`grep '^e[[:space:]]' $temp_file | wc | awk '{print $1}'`
            edge_string=`echo $desired_num_edges | awk '{printf "%06d",$1}'`
            out_file=${graph_type}_${node_string}_${edge_string}_${weight_string}_${seed_string}.gph
            echo "c $num_nodes $true_num_edges $max_weight $seed" > $out_file
            sed "s/^g[[:space:]][0-9]*[[:space:]][0-9]*/g $num_nodes $true_num_edges/" $temp_file >> $out_file
            rm $temp_file
        fi
        echo "successfully created graph in $out_file"
        successful_attempts=$((successful_attempts + 1))
    fi
    seed=$((seed + 1))
done

echo "### Created $num_instances instances in $attempts attempts ###"
