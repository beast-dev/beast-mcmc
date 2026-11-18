#!/bin/bash
#
# sign_jar_natives.sh <signing identity> <directory or jar>...
#
# Signs the macOS native libraries (*.dylib, *.jnilib) that are packed inside jar files, in place.
#
# Apple's notary service looks inside jars and rejects a package if any Mach-O binary in one is not
# signed with a Developer ID certificate, the hardened runtime and a secure timestamp. Jars that carry
# their own native code (e.g. zstd-jni, which is also merged into beast.jar) arrive with ad hoc or no
# signatures, so they have to be signed here before the .app bundles and the DMG are signed.
#
# Set CODESIGN_TIMESTAMP=--timestamp=none to sign without a timestamp server (for testing only).

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "usage: $0 <signing identity> <directory or jar>..." >&2
    exit 1
fi

SIGNER="$1"
shift
TIMESTAMP="${CODESIGN_TIMESTAMP:---timestamp}"

while IFS= read -r JAR; do
    ENTRIES=$(unzip -Z1 "$JAR" | grep -E '\.(dylib|jnilib)$' || true)
    if [ -z "$ENTRIES" ]; then
        continue
    fi

    JAR_PATH="$(cd "$(dirname "$JAR")" && pwd)/$(basename "$JAR")"
    WORK=$(mktemp -d)

    while IFS= read -r ENTRY; do
        unzip -q -o "$JAR_PATH" "$ENTRY" -d "$WORK"
        codesign --force --options runtime "$TIMESTAMP" --sign "$SIGNER" "$WORK/$ENTRY"
        echo "signed $ENTRY in $JAR"
    done <<< "$ENTRIES"

    # put the signed libraries back in the jar, replacing the originals
    (cd "$WORK" && zip -q "$JAR_PATH" -@ <<< "$ENTRIES")
    rm -rf "$WORK"
done < <(find "$@" -type f -name "*.jar")
