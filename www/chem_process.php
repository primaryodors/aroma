<?php

function chem_python(): string
{
    $root = dirname(__DIR__);
    $configured = getenv("AROMA_PYTHON");
    if ($configured) $candidates = [$configured];
    elseif (is_file("$root/venvname") && ($environment = trim(file_get_contents("$root/venvname"))) !== "")
    {
        if ($environment[0] !== "/") $environment = "$root/$environment";
        $candidates = ["$environment/bin/python3", "$environment/bin/python"];
    }
    else
    {
        $candidates = ["/usr/local/bin/python3", "/usr/bin/python3"];
        foreach (explode(PATH_SEPARATOR, getenv("PATH") ?: "") as $directory)
            if ($directory !== "") $candidates[] = "$directory/python3";
    }
    foreach ($candidates as $candidate)
        if (is_file($candidate) && is_executable($candidate)) return $candidate;
    throw new RuntimeException("The configured Python interpreter is unavailable.", 503);
}

function chem_process_result(array $arguments, ?string $directory = null, float $timeout = 20.0): array
{
    foreach (["proc_open", "proc_get_status", "proc_terminate", "proc_close"] as $function)
        if (!function_exists($function))
            throw new RuntimeException("Process execution is disabled on this host ($function).", 503);
    $command = PHP_VERSION_ID >= 70400 ? $arguments
        : "exec ".implode(" ", array_map("escapeshellarg", $arguments));
    $process = @proc_open($command, [
        0 => ["file", "/dev/null", "r"],
        1 => ["pipe", "w"],
        2 => ["redirect", 1],
    ], $pipes, $directory);
    if (!is_resource($process))
        throw new RuntimeException("The host could not start the molecular tool.", 503);
    stream_set_blocking($pipes[1], false);
    $output = "";
    $began = microtime(true);
    $status = -1;
    try
    {
        while (true)
        {
            $output .= stream_get_contents($pipes[1]);
            if (strlen($output) > 262144)
                throw new RuntimeException("Molecular tool output exceeded its limit.", 422);
            $state = proc_get_status($process);
            if (!$state["running"])
            {
                $status = $state["exitcode"];
                $output .= stream_get_contents($pipes[1]);
                break;
            }
            if (microtime(true)-$began > $timeout)
                throw new RuntimeException("Molecule generation timed out.", 504);
            usleep(20000);
        }
    }
    finally
    {
        if (isset($state) && $state["running"]) proc_terminate($process, 9);
        fclose($pipes[1]);
        $closed = proc_close($process);
        if ($status < 0) $status = $closed;
    }
    return ["status" => $status, "output" => $output];
}

function run_chem_process(array $arguments): int
{
    try
    {
        $result = chem_process_result($arguments);
        echo $result["output"];
        if ($result["status"] !== 0) http_response_code(422);
        return $result["status"];
    }
    catch (Throwable $error)
    {
        $code = $error->getCode();
        http_response_code(in_array($code, [422, 503, 504], true) ? $code : 500);
        echo $error->getMessage();
        error_log("SMILES conversion: ".$error->getMessage());
        return -1;
    }
}
