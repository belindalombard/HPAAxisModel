function check_parameters(required_params)
    for i = 1:length(required_params)
        param = required_params{i};
        if ~evalin('base', sprintf('exist(''%s'', ''var'')', param))
            error('The variable "%s" does not exist. Please run the entire pipeline to ensure that the variables are created in the scope.', param);
        end
    end
end