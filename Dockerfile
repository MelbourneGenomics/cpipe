FROM ubuntu

# Start using bash
SHELL ["/bin/bash", "-c"]

# Install basic linux tools that we need to make python
RUN apt-get update && apt-get install -y curl make build-essential libsqlite3-dev

# Copy in the repository
ADD . /opt/cpipe

# Move into the cpipe dir
WORKDIR /opt/cpipe

# Run the install script
RUN ./install.sh

# Run the main script
ENTRYPOINT ["./cpipe"]
