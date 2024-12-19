# About

This code is a port of a port of a port, but it seems to work. If anyone
that actually understands the math wants to takeover this codebase, please let me know.

The Deriche filter is an order-K approximation of the Gaussian kernel, which is
combined with linear binning to produce a fast and accurate implementation of KDE.
The naive KDE algorithm is O(n^2), while this implementation is O(n+m) wherem is the
number of bins.


# References

This code is a port of this Python implementation of the Deriche approximation KDE algorithm found here:
https://github.com/liuzh-buaa/fast_kde

Which is itself a port of Javascript implementation... which was the implementation of a paper
that was based on a port of the Deriche computer vision kernel implemented in C by Getreuer.

The original Javascript implementation is here:
https://github.com/uwdata/fast-kde

The paper on the implementation for KDE can be found here:
https://idl.uw.edu/papers/fast-kde
