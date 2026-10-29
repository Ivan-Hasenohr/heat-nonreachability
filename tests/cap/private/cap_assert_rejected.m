function cap_assert_rejected(operation,message_fragment)
% The validation error must mention the intended input contract.
rejected = false;
try
    operation();
catch failure
    rejected = contains(failure.message,message_fragment);
end
assert(rejected,'Invalid input was not rejected by the expected validation.');
end
