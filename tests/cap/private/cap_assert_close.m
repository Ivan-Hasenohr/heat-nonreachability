function cap_assert_close(actual,expected,tolerance)
% Compare narrow, independently computed interval references.
% Overlap is a regression check, not a proof of full interval containment.
if nargin < 3, tolerance = 1e-9; end
actual = intval(actual); expected = intval(expected);
assert(isequal(size(actual),size(expected)),'Unexpected output shape.');
assert(all(isfinite(inf(actual(:)))) && all(isfinite(sup(actual(:)))));
assert(all(inf(actual(:)) <= sup(expected(:))) ...
    && all(sup(actual(:)) >= inf(expected(:))), ...
    'Computed interval and independent reference are disjoint.');
assert(all(sup(abs(actual(:)-expected(:))) < tolerance), ...
    'Unexpected discrepancy or interval width on this small case.');
end
