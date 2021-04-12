
# Run tests with coverage
run(`julia --code-coverage=user --inline=yes --project -e 'cd("./test"); include("runtests.jl")'`)

using Coverage
# process '*.cov' files
coverage = process_folder(); # defaults to src/; alternatively, supply the folder name as argument
# Get total coverage for all Julia files
covered_lines, total_lines = get_summary(coverage);

covered_lines/total_lines*100

# Clean up!!
clean_folder("src")
clean_folder("test")