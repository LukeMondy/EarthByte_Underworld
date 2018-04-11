
imgname := olympusmonds/earthbyte_underworld
cntname := ebuw


.PHONY: build
build:
	docker build \
		-t $(imgname) \
		-f ebuw.dockerfile \
		.

.PHONY: run
run: 
	docker run \
		-u "$(shell id -u)" \
		-w "$(shell echo $$HOME)" \
		--hostname EarthByteUnderworld \
		-it \
		--rm \
		-v /etc/passwd:/etc/passwd:ro \
		-v /etc/group:/etc/group:ro \
		-v /etc/localtime:/etc/localtime:ro \
		-v /home:/home \
		-v /:/host \
		$(imgname) bash
