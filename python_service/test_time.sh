start=$(date +%s)

# 1. Submit job and capture job_id
job_id=$(curl -X POST http://localhost:8000/v1/sdf/descriptors \
  -F "file=@/home/kapnevets/Projects/rdkitwebserver/python_service/examples/Compound_000000001_000500000.sdf" | jq -r '.job_id')

# 2. Poll until job is finished
while true; do
  status=$(curl -s http://localhost:8000/v1/job/$job_id | jq -r '.status')
  if [ "$status" == "done" ]; then
    break
  fi
  sleep 1
done

end=$(date +%s)
echo "Total elapsed: $((end - start)) seconds"