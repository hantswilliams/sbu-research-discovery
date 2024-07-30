## Building Docker Image

To build the docker image, run the following command:
- Image-name here should be `shp-researcher`
- The buildx linx/amd64 is used to build the image for linux/amd64 platform which is what we mostly have/use in the cloud.

```bash
docker buildx build --platform linux/amd64 -t <image-name> .
```

So it could be:
```bash
docker buildx build --platform linux/amd64 -t shp-researcher .
```

## Running Docker Image

To run the docker image, run the following command:

```bash
docker run -p 5009:5009 <image-name>
```

So it could be:
```bash
docker run -p 5009:5009 shp-researcher
```

## Tagging Docker Image

To tag the docker image, run the following command:

```bash
docker tag <image-name> <docker-hub-username>/<image-name>
```

So it could be:
```bash
docker tag shp-researcher hants/shp-researcher
```

## Pushing Docker Image

To push the docker image to docker hub, run the following command:

```bash
docker push <docker-hub-username>/<image-name>
```

So it could be:
```bash
docker push hants/shp-researcher
```