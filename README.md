# Planar-Elastic-Metrics

This is a numerical implementation of the work in 

@article{kurtek2018simplifying,
  title={Simplifying transforms for general elastic metrics on the space of plane curves},
  author={Kurtek, Sebastian and Needham, Tom},
  journal={arXiv preprint arXiv:1803.10894},
  year={2018}
}

The transforms are normalized to have a=1, so that the 1-parameter family of metrics considered are the elastic metrics g^{1,b}. Due to non-injectivity and numerical issues outlined in the paper, there are two regimes where the implementation has different performance:

For “large” b values, the transform works in a straightforward and fast way. “Large” depends on how complicated the curves are. Usually b>0.4 works, lower values work for very simple curves, and b>0.5 always works for theoretical reasons.

For “small” b values, we use the transform plus a path straightening algorithm to fix the noninjectivity issues. The algorithm is simplistic, but slow in its current form.

To run the program, first compile DynamicProgrammingQ.c as a mex file. Then use main_large_b(c1,c2,b) or main_small_b(c1,c2,b). Typing “example” will automatically run an example with a large b value. Example curves are provided in the comp_curves.mat file.

Please feel free to contact me with any questions or to report any issues.
