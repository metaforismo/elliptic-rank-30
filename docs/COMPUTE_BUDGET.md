# Compute allocation policy

Expensive computation is released in stages:

1. exact cheap algebra and duplicate rejection;
2. local and finite-field filters;
3. descent and visible-point searches;
4. incremental independence checks;
5. height and saturation only for promoted candidates.

Every batch reports cost per new independent point and the dominant failure
mode.  Bounds are increased only after the failure analysis justifies it.
