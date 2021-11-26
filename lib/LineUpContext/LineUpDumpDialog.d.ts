import { ConnectedProps } from "react-redux";
declare const connector: import("react-redux").InferableComponentEnhancerWithProps<{
    setLineUp_dump: (dump: any) => any;
    setLineUp_visibility: (vis: any) => any;
}, {}>;
declare type PropsFromRedux = ConnectedProps<typeof connector>;
declare type Props = PropsFromRedux & {
    openDialog: any;
    setOpenDumpDialog: any;
};
export declare const LineUpDumpDialog: import("react-redux").ConnectedComponent<({ openDialog, setOpenDumpDialog, setLineUp_dump, setLineUp_visibility, }: Props) => JSX.Element, import("react-redux").Omit<{
    setLineUp_dump: (dump: any) => any;
    setLineUp_visibility: (vis: any) => any;
} & {
    openDialog: any;
    setOpenDumpDialog: any;
}, "setLineUp_dump" | "setLineUp_visibility">>;
export {};
//# sourceMappingURL=LineUpDumpDialog.d.ts.map