## How to use BEAST checkpointing

To write a state file to disk at a given iteration number, add the following program argument: **'-save_at iteration#'** (without the quotation marks).
This will create a time-stamped file 'beast_state_xxxx.xx.xx' containing all the necessary information.

To write a state file to disk every x iterations, add the following program argument: **'-save_every x'** (without the quotation marks).
If no file name has been provided and there is at least 1 second time difference between the previous attempt to write to file, a new time-stamped file will be created every time the information needs to be written to disk.

To write the state file to a specific file name, add the following program argument: **'-save_state filename'** (without the quotation marks).
Every time a new state needs to be stored, the file's content will be overwritten.

To load a BEAST state file from disk and resume a previous analysis, add the following program argument: **-'load_state filename'** (without the quotation marks).

