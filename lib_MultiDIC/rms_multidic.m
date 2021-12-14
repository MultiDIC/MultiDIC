function output = rms_multidic(input)
        output = sqrt(mean(input .* input));
end
