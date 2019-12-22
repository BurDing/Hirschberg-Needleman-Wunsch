for ((i = 0; i < 10; i += 1))
do
  for ((len = 0; len <= 9000; len += 500))
  do
    echo $i, $len
    python3 memory.py --method global --length $len
    python3 memory.py --method hirschberg --length $len
  done
done
