
This is a "to-do" list for the future. 

* remove the step cloning the `dismat` matrix for use with the customized weight function; it should be possible to use a reference to the original dismat instead, without copying.
* right now we use iterators to find "good" simplices to index the rows and columns of our matrices.  these iterators would be more efficient if we skipped simplices that are too big.  we may need to change something about the way the iterators are structured in general in order to do this, however.  it may not be possible to implement this in optimal reps until this more fundamental change is made to the exhact library.
