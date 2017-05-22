## How to use BEAST debugging

To write a dump file to disk at a given iteration number, add the following program argument: **'-dump_state iteration#'** (without the quotation marks).
This will create a time-stamped file 'beast_debug_xxxx.xx.xx' containing all the necessary information.

To write a dump file to disk every x iterations, add the following program argument: **'-dump_every x'** (without the quotation marks).
If no file name has been provided and there is at least 1 second time difference between the previous attempt to write to file, a new time-stamped file will be created every time the information needs to be written to disk.

To write the dump file to a specific file name, add the following program argument: **'-save_dump filename'** (without the quotation marks).
Every time a new state needs to be stored, the file's content will be overwritten.

To load a BEAST debugging file from disk and resume a previous analysis, add the following program argument: **-'load_dump filename'** (without the quotation marks).

