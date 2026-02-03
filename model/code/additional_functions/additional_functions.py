def merge_parameters_stressor(arg_params, json_params, default_params):
    merged_params = {
        key: arg_params.get(key) or json_params.get(key) or default_params[key]
        for key in default_params
    }

    return merged_params