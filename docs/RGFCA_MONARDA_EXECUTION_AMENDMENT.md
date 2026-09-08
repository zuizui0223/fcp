# Monarda: pre-outcome execution clarification

The original contract at `336b4b8` is preserved byte-for-byte in Git. Before
any Monarda pixel opening, both execution tasks confirmed no JPEG decoding,
polygon rasterization or model loading. The peer explicitly handed Monarda
editing and execution to this desktop task; its sharedness-v2 work is separate.

## What is corrected

The original contract specifies `FrozenFlowerColourEstimator.measure()` but
also prohibits colour measurement. That exact method computes CIELAB internally
for original/flip measurement and admission. Therefore, “never measures colour”
was inaccurate. The [execution amendment](supporting/rgfca_monarda_execution_amendment_v1.json)
explicitly permits this incidental computation, retaining the originally
specified full method and failure semantics. It does not save or analyze
continuous colour, filter the scored images by colour admission, or validate
colour accuracy. No alternative mask algorithm is introduced.

The 110 images, 109 annotated-image scoring set, unknown image, model weights,
runtime, rasterization and `.70/.35/.70` operational floors are unchanged.
JRC object/box-derived floors are not independently established petal-pixel
accuracy criteria. Success remains limited generic flower-region agreement.

## Before images

The execution guard must verify a separately committed authorization referring
to qualified code and successful CI, exact code/contract/model identities,
the pinned archive/member census and the prospectively specified Windows
environment. The isolated immutable runtime cache does not modify the repo's
current or historical runtime. Text-only identity comparisons explicitly
normalize CRLF to LF; model/archive identities remain exact byte hashes.

Only the canonical output directory may be used. It must not exist. Its durable
start receipt and complete 110-row pending ledger are created before model
setup or pixels. Alignment and scoring are journalled per image, with atomic
complete-ledger checkpoints. A failure keeps incomplete/unknown rows and cannot
produce a successful-subset conclusion. Existing output blocks every rerun;
recovery requires a separately recorded decision, not silent overwriting.

This document is not an execution authorization and does not report image
performance. The next gate is implementation qualification and a separate
authorization bound to the exact qualified commit. Ecological results remain
unchanged, including the reserve flower-specific differential `p = 0.087`.

## Qualification and single-execution authorization

The amended implementation at `1b6ccc7987baaa02edf8fea8c5af1f0f09c144f7`
passed [qualification run 34193780128](https://github.com/zuizui0223/fcp/actions/runs/34193780128):
23 focused artificial tests, exact runtime/model identity checks and the
pre-outcome firewall. Local qualification including intake tests passed 53/53.
The [separate authorization](supporting/rgfca_monarda_execution_authorization_v1.json)
binds the exact qualified source blobs and permits the one guarded run after
local environment/input checks. No real image was opened to choose this
implementation, environment or gate. Actual performance is not reported here.
