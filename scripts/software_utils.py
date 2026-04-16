def get_support_software_list(event_type):
    """
    Return the list of supported software for a given splicing event type.
    """

    if event_type in ("SE", "A3SS", "A5SS", "RI"):
        return ["SUPPA2", "rMATS", "PSI-Sigma", "MAJIQ", "Spladder", "Whippet"]

    elif event_type in ("AF", "AL"):
        return ["SUPPA2", "MAJIQ", "Whippet"]

    elif event_type in ("MX",):
        return ["SUPPA2", "rMATS", "PSI-Sigma", "MAJIQ", "Spladder"]

    else:
        return []
