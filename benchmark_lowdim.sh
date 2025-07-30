echo "Generating initial points"
python3 ./scripts/lowdim_gen_init.py
echo "Benchmarking scipy optim"
python3 ./scripts/lowdim_bmrk.py > outputs/py_bmrk_table.md
echo "Benchmarking fdaPDE optim"
./build/advanced_optim_applic output_csv lowdim
echo "Generating plots"
python3 ./scripts/lowdim_plot.py