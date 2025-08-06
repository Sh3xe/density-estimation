echo "Generating initial points"
python3 ./scripts/lowdim_gen_init.py
echo "Benchmarking scipy optim"
python3 ./scripts/lowdim_bmrk.py > outputs/py_bmrk_table.md
echo "Benchmarking fdaPDE optim"
./build/advanced_optim_applic lowdim output_csv
echo "Generating plots"
python3 ./scripts/lowdim_plot.py