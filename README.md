# SubgraphCount
To build the project run in the same directory 


$cmake CMakeLists.txt

$make 

It will create the following executable:
1. SubgraphCount.out
2. FormatConverter.out
3. sanitize.out
4. ExactCount.out
5. GraphProperties.out
6. BaselineComparison.out

To run TETRIS, we use SubgraphCount.out.

## Input file format:

SubgraphCount.out expects an input file in the CSR format. 
If the input file in the raw edges format, 
the following steps will convert it into the CSR format.

#### Convert input files with raw edge data into .edge format
First run the sanitize.out on the input file
to clean the raw data. 
It takes as input a file with edges list as 
pairs of strings, and creates a file in the .edges format. 

Usage: ./sanitize.out path/to/out/directory/ input_filename

The .edges format has the first line with the number 
of nodes and edges. Each line has a distinct edge with 
node labels as ints starting from 0.

#### Convert .edges file into the CSR format

The second step is to take this file in the .edges
format and convert it into a .CSR format. The
GraphFormatConverter.out executable does this work.

Usage: ./FormatConverter.out path/to/out/directory/ input_filename

Input file must be in .edge format, 
where first line has number of nodes and edges, 
and every line has a distinct undirected edge.


## Running TETRIS

To run our algorithm TETRIS, run SubgraphCount.out.

Usage: ./SubgraphCount script_file

The script_file contains all the input configuration parameters. 
See sample_config.txt for parameter details. 
