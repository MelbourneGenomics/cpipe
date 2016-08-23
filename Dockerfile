FROM ubuntu
SHELL ["/bin/bash"]
ADD .
RUN install/dependencies/ubuntu_16.04.sh
RUN install/download_assets.sh
RUN install/install.sh
ENTRYPOINT docker/entrypoint.sh
