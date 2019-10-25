
imgname := lmondy/earthbyte_underworld

WORKDIR ?= $(shell pwd)

.PHONY: build
build:
	docker build \
		-t $(imgname) \
		-f ../docker-base/Dockerfile \
		.

.PHONY: run
run: 
	docker run \
		-u "$(shell id -u)" \
		-w /models \
		--hostname EarthByteUnderworld \
		-it \
		--rm \
		-v /etc/passwd:/etc/passwd:ro \
		-v /etc/group:/etc/group:ro \
		-v /etc/localtime:/etc/localtime:ro \
		-v /home:/home \
		-v $(WORKDIR):/host \
		$(imgname) bash
