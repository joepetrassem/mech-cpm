"""Turn a MATLAB JUnit report into GitHub Actions annotations.

Each failed or errored test becomes an error annotation (visible on the
pull request and in the run summary) with its diagnostic message, and a
notice gives the totals.
"""
import sys
import xml.etree.ElementTree as ET


def escape(text):
    # Workflow-command escaping: % first, then newlines.
    return text.replace("%", "%25").replace("\r", "").replace("\n", "%0A")


def main(path):
    try:
        root = ET.parse(path).getroot()
    except (OSError, ET.ParseError) as exc:
        print(f"::warning title=Test report::Could not read {path}: {exc}")
        return
    cases = root.iter("testcase")
    total = failed = 0
    for case in cases:
        total += 1
        problem = case.find("failure")
        if problem is None:
            problem = case.find("error")
        if problem is None:
            continue
        failed += 1
        name = f"{case.get('classname', '')}/{case.get('name', '')}"
        message = (problem.get("message") or "") + "\n" + (problem.text or "")
        print(f"::error title={escape(name)}::{escape(message.strip()[:3000])}")
    print(f"::notice title=Test summary::{total} tests, {failed} failed or errored")


if __name__ == "__main__":
    main(sys.argv[1])
