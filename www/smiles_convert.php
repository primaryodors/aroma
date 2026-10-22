<?php
require_once(__DIR__."/chem_process.php");

function viewer_tool(string $name, string $setting = ""): string
{
    $root = dirname(__DIR__);
    $configured = $setting ? getenv($setting) : false;
    $candidates = $configured ? [$configured] : [
        "$root/$name", "/usr/local/bin/$name", "/usr/bin/$name"
    ];
    if (!$configured)
        foreach (explode(PATH_SEPARATOR, getenv("PATH") ?: "") as $directory)
            if ($directory !== "") $candidates[] = "$directory/$name";
    foreach ($candidates as $candidate)
        if (@is_file($candidate) && @is_executable($candidate)) return $candidate;
    throw new RuntimeException("The host cannot find or execute $name. Check its installation and executable path.", 503);
}

function viewer_sdf(string $path): string
{
    if (!is_file($path) || filesize($path) > 2097152)
        throw new RuntimeException("The molecular tool did not produce a valid structure.", 422);
    $sdf = file_get_contents($path);
    $lines = preg_split('/\r?\n/', $sdf);
    $atoms = isset($lines[3]) ? (int)substr($lines[3], 0, 3) : 0;
    if ($atoms < 1 || $atoms > 512 || count($lines) < $atoms+4 || strpos($sdf, "M  END") === false)
        throw new RuntimeException("The molecular tool returned an empty or unsupported structure.", 422);
    $positions = [];
    for ($i=4; $i<$atoms+4; $i++) $positions[substr($lines[$i], 0, 30)] = true;
    if ($atoms > 1 && count($positions) < 2)
        throw new RuntimeException("The molecular tool returned a structure without coordinates.", 422);
    return $sdf;
}

function viewer_python(): string
{
    if (getenv("AROMA_PYTHON")) return viewer_tool("python3", "AROMA_PYTHON");
    $root = dirname(__DIR__);
    if (is_file("$root/venvname"))
    {
        $environment = trim(file_get_contents("$root/venvname"));
        if ($environment !== "")
        {
            if ($environment[0] !== "/") $environment = "$root/$environment";
            foreach (["python3", "python"] as $name)
                if (is_executable("$environment/bin/$name")) return "$environment/bin/$name";
            throw new RuntimeException("The configured Python environment has no executable interpreter.", 503);
        }
    }
    return viewer_tool("python3");
}

function convert_viewer_molecule(array $request): array
{
    $smiles = $request["smiles"] ?? null;
    if (!is_string($smiles) || trim($smiles) === "" || strlen($smiles) > 4096 || strpos($smiles, "\0") !== false)
        throw new RuntimeException("Enter a valid SMILES string.", 400);
    if (strpos($smiles, "{") !== false || strpos($smiles, "}") !== false)
    {
        $plain = preg_replace('/\{[A-Za-z][A-Za-z0-9_]{0,13}\}/', '', $smiles);
        if (strlen($smiles) > 255 || strpbrk($plain, "{}") !== false)
            throw new RuntimeException("Custom atom labels accept up to 14 letters or digits in a short molecular input.", 400);
    }
    $transforms = [];
    foreach (["brot", "rflp"] as $option)
    {
        $value = $request[$option] ?? "";
        if (!is_string($value) || strlen($value) > 4096 ||
            ($value !== "" && !preg_match('/\A[A-Za-z0-9_.+\s-]+\z/', $value)))
            throw new RuntimeException("Invalid transformation arguments.", 400);
        $transforms[$option] = preg_split('/\s+/', trim($value), -1, PREG_SPLIT_NO_EMPTY);
    }
    $rotation = $transforms["brot"];
    if ($rotation && $rotation[0] === "minc") array_shift($rotation);
    if (count($rotation) % 3 !== 0)
        throw new RuntimeException("Bond rotations accept atom, atom, angle triples.", 400);
    for ($i=2; $i<count($rotation); $i+=3)
        if (!is_numeric($rotation[$i]) || !is_finite((float)$rotation[$i]))
            throw new RuntimeException("A rotation angle is invalid.", 400);

    $root = dirname(__DIR__);
    $temporary = getenv("AROMA_TMPDIR") ?: sys_get_temp_dir();
    $workspace = @tempnam($temporary, "aroma-view-");
    if ($workspace === false)
        throw new RuntimeException("The host has no writable temporary directory.", 503);
    unlink($workspace);
    if (!@mkdir($workspace, 0700))
        throw new RuntimeException("The host could not create a temporary workspace.", 503);
    $output = "$workspace/molecule.sdf";
    $log = "";
    $run = function(array $command) use ($root, &$log) {
        $result = chem_process_result($command, $root);
        $log .= $result["output"];
        if ($result["status"] !== 0)
            throw new RuntimeException("Molecular tool failed: ".trim($result["output"]), 422);
    };
    try
    {
        if (strpos($smiles, "{") !== false)
            $run([viewer_tool("test/mol_assem_test"), $smiles, $output]);
        else
            $run([viewer_python(), "$root/data/viewer_sdf.py", $smiles, $output]);
        viewer_sdf($output);
        if ($transforms["brot"])
        {
            $rotated = "$workspace/rotated.sdf";
            $run(array_merge([viewer_tool("test/bond_rotation_test"), $output],
                             $transforms["brot"], ["--output", $rotated]));
            viewer_sdf($rotated);
            if (!rename($rotated, $output)) throw new RuntimeException("Could not store the rotated structure.", 500);
        }
        if ($transforms["rflp"])
            $run(array_merge([viewer_tool("bin/ringflip"), $output], $transforms["rflp"]));
        return ["sdf" => viewer_sdf($output), "log" => trim(str_replace($workspace, "[temporary files]", $log))];
    }
    finally
    {
        foreach (glob("$workspace/*") ?: [] as $file) if (is_file($file)) unlink($file);
        rmdir($workspace);
    }
}
