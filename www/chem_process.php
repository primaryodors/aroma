<?php

function run_chem_process(array $arguments): int
{
    $process = proc_open($arguments, [
        0 => ["file", "/dev/null", "r"],
        1 => ["pipe", "w"],
        2 => ["redirect", 1],
    ], $pipes);
    if (!is_resource($process))
    {
        http_response_code(500);
        return -1;
    }
    fpassthru($pipes[1]);
    fclose($pipes[1]);
    $status = proc_close($process);
    if ($status !== 0) http_response_code(422);
    return $status;
}
