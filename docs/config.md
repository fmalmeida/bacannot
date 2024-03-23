# Configuration File

To download a configuration file template users just use `--get_config` parameter. Using a config file your code is lot more clean and concise.

```bash
# get config
nextflow run fmalmeida/bacannot --get_config
# run with config
nextflow run fmalmeida/bacannot -c [path-to-config]
```

Default configuration
---------------------

```groovy
{% include 'defaults.config' %}
```