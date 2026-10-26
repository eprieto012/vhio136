ENV_REPRO=/mnt/bioinfnas/immuno/eprieto/Rmagic/miniforge3/envs/renv_vhio136
PROJECT=/mnt/bioinfnas/immuno/eprieto/analyses/vhio136/renv_vhio136

export PATH="$ENV_REPRO/bin:$PATH"
export PKG_CONFIG_LIBDIR="$ENV_REPRO/lib/pkgconfig:$ENV_REPRO/share/pkgconfig"
export PKG_CONFIG_PATH="$PKG_CONFIG_LIBDIR"

cd "$PROJECT"

env -u R_LIBS -u R_LIBS_SITE -u R_LIBS_USER "$ENV_REPRO/bin/R"