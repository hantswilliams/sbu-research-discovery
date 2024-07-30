FROM python:3.11.6-alpine3.18

## install necessary required for gcc
RUN apk add --no-cache gcc musl-dev linux-headers

WORKDIR /app
COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt
COPY . .
EXPOSE 5009
CMD [ "python", "app.py" ]